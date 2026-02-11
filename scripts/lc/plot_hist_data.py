"""
This file gathers the historical timings for certain tests and creates a gitlab pages
HTML file. This file should be run from a virtual python environment, which can be
created using the following commands:

python3 -m venv --system-site-packages /path/to/env
source /path/to/env/bin/activate
# or if using csh
# source /path/to/env/bin/activate.csh
python3 -m pip install --progress-bar off -U pip
python3 -m pip install --progress-bar off caliper-reader pyyaml pandas plotly mpi4py
An environment should already exist in /g/g15/sphapp/docsenv
"""

import os, sys, yaml, glob, time, shutil, argparse
from datetime import datetime
import caliperreader as cr
import multiprocessing as mp
import concurrent.futures
import pandas as pd
import numpy as np
import plotly.express as px
from mpi4py import MPI
import hist_data_utils as hdu

# How many months to plot
num_of_months = 6
# Which region to plot
plot_region = "advance"
colors = px.colors.qualitative.Set1

def init_worker():
    global reg_names, spot_metric
    reg_names = ["advance", "ConnectivityMap_computeConnectivity", "computeVoronoiVolume"]
    spot_metric = "avg#inclusive#sum#time.duration"

def extract_data(cali_file):
    global reg_names, spot_metric
    (records, gls) = cr.read_caliper_contents(cali_file)
    runday = int(gls["launchdate"])
    rtimes = {}
    for rec in records:
        if ('path' in rec):
            for rg in reg_names:
                if (rec["path"][-1] == rg):
                    tval = float(rec[spot_metric])
                    if (rg in rtimes.keys()):
                        rtimes[rg] += tval
                    else:
                        rtimes[rg] = tval
    return runday, rtimes

def get_spec(cali_file):
    # If cali_file is just a dir, grab a cali file from latest
    if os.path.isdir(cali_file):
        cali_file = glob.glob(os.path.join(cali_file, "latest/*.cali"))[0]
    gls = cr.read_caliper_globals(cali_file)
    return gls["spec"]

def compare_time(cfile, num_of_months):
    tdelta = (time.time() - os.path.getmtime(cfile))/ (60*60*24)
    if (tdelta < 30*num_of_months):
        return True
    return False

def get_latest_files(cdir, test_name, num_of_months):
    all_files = glob.glob(os.path.join(cdir, f"*/{test_name}*.cali"))
    cali_files = [cf for cf in all_files if compare_time(cf, num_of_months)]
    return cali_files

# Convert array of dictionaries into panda dataframe
def convert_to_dataframe(in_array):
    data = []
    for in_data in in_array:
        mac_name = in_data["machine"]
        spec_name = in_data["spec"]
        config = in_data["config"]
        data_dict = in_data["data_dict"]
        for test_name, tdata in data_dict.items():
            dates = tdata["dates"]
            times = tdata["times"]
            for date, time_dict in zip(dates, times):
                for reg, val in time_dict.items():
                    data.append({"Machine": mac_name,
                                 "Spec": spec_name,
                                 "Config": config,
                                 "Test Name": test_name,
                                 "Region": reg,
                                 "Date": str(date),
                                 "Time": val})
    df = pd.DataFrame(data)
    df = df.sort_values(by="Date")
    min_per_day = df.groupby(["Config", "Test Name", "Date"])["Time"].transform("min")
    max_per_day = df.groupby(["Config", "Test Name", "Date"])["Time"].transform("max")
    df["extreme"] = "other"
    df.loc[df["Time"].eq(min_per_day), "extreme"] = "Min"
    df.loc[df["Time"].eq(max_per_day), "extreme"] = "Max"
    return df

def create_plot(in_df):
    mac_name = in_df["Machine"].unique()[0]
    specs = in_df["Spec"].unique()
    test_names = in_df["Test Name"].unique()
    num_tests = len(test_names)
    figs = []
    # Pick the region to plot
    df = in_df[in_df["Region"] == plot_region]
    for test in test_names:
        test_df = df[df["Test Name"] == test]
        plot_df = test_df[test_df["extreme"].isin(["Min"])]
        fig = px.line(title=test)
        for cindx, spec in enumerate(specs):
            ccolor = colors[cindx]
            cdf = test_df[test_df["Spec"] == spec]
            # Plot the min, max and mean for any runs on a given day
            mean_df = cdf.groupby("Date")["Time"].mean().reset_index()
            min_df = cdf.groupby("Date")["Time"].min().reset_index()
            max_df = cdf.groupby("Date")["Time"].max().reset_index()
            fig.add_scatter(x=max_df["Date"],
                            y=max_df["Time"],
                            legendgroup=spec,
                            showlegend=False,
                            line=dict(color=ccolor))
            fig.add_scatter(x=min_df["Date"],
                            y=min_df["Time"],
                            legendgroup=spec,
                            showlegend=False,
                            fill="tonexty",
                            line=dict(color=ccolor))
            fig.add_scatter(x=mean_df["Date"],
                            y=mean_df["Time"],
                            line=dict(color=ccolor),
                            marker=dict(color=ccolor, size=10),
                            mode="lines+markers",
                            name=f"{spec}",
                            legendgroup=spec,
                            showlegend=True)

        figs.append(fig)
    return figs

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument("--doc-dir", default="public")
    parser.add_argument("--bench", default="/usr/workspace/sduser/Spheral/benchmarks")
    parser.add_argument("--threads", default=None, type=int)
    parser.add_argument("--pkg-name", type=str, help="Either Spheral or LLNLSpheral")
    parser.add_argument("--test", action="store_true")
    args = parser.parse_args()
    cur_dir = os.getcwd()
    doc_dir = os.path.join(cur_dir, args.doc_dir)
    if (args.pkg_name == "Spheral"):
        sys.append("../../tests/performance")
        import perf_tests as pt
    else:
        sys.append("../../../tests/performance")
        import llnl_perf_tests as pt
    if (rank == 0):
        if (os.path.exists(doc_dir)):
            shutil.rmtree(doc_dir)
        os.makedirs(doc_dir)
    test_names = pt.get_all_test_names()
    if (args.test):
        test_names = [test_names[0]]
    # Benchmark directory is organized as /configs/machines/run_dates/*.cali
    # Each set of caliper files contain a single spec
    # Multiple configs could span a single machine+spec
    config_dirs = glob.glob(os.path.join(args.bench, "*"))
    # Extract all machine directories
    all_mac_specs = {}
    all_specs = {}
    if (rank == 0):
        for cc in config_dirs:
            mac_dirs = glob.glob(os.path.join(cc, "*"))
            for mc_path in mac_dirs:
                mac_name = os.path.basename(mc_path)
                cspec = get_spec(mc_path)
                key = mac_name+" "+cspec
                all_specs[key] = cspec
                if key in all_mac_specs:
                    all_mac_specs[key].append(mc_path)
                else:
                    all_mac_specs.update({key: [mc_path]})
    all_mac_specs = comm.bcast(all_mac_specs, root=0)
    all_specs = comm.bcast(all_specs, root=0)
    all_data = []
    # Loop over each machine+spec and extract the data
    for mac_spec, mc_dirs in all_mac_specs.items():
        mac_name = os.path.basename(mc_dirs[0])
        spec_name = all_specs[mac_spec]
        config_str = mac_spec.replace(" ", "_")
        timer_data = {}
        my_tests = test_names[rank::size]
        for test_name in my_tests:
            test_files = []
            print(f"Rank {rank} processing {mac_spec} {test_name} Caliper files")
            for mc in mc_dirs:
                test_files.extend(get_latest_files(mc, test_name, num_of_months))
            with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads, initializer=init_worker) as pool:
                result = list(pool.map(extract_data, test_files))
            ctftimes = []
            xaxis = []
            for i in result:
                xaxis.append(datetime.fromtimestamp(i[0]).replace(hour=0, minute=0, second=0, microsecond=0))
                ctftimes.append(i[1])
            timer_data.update({test_name: {"dates": xaxis, "times": ctftimes}})
        gtimer_data = comm.gather(timer_data, root=0)
        if (rank == 0):
            # Flatten list tests
            lt_data = {k: v for d in gtimer_data for k, v in d.items()}
            all_data.append({"config": config_str, "machine": mac_name,  "spec": spec_name, "data_dict": lt_data})
        comm.Barrier()
    if (rank == 0):
        df_data = convert_to_dataframe(all_data)
        macs = df_data["Machine"].unique()
        file_link_dict = {mac: f"{mac}.html" for mac in macs}
        file_link_dict.update({"Home": "index.html"})
        for mac in macs:
            mac_df = df_data[df_data["Machine"] == mac]
            figs = create_plot(mac_df)
            page_content = hdu.create_page_content(figs)
            hdu.create_html_file(os.path.join(doc_dir, f"{mac}.html"), file_link_dict, mac, page_content)
        hdu.create_html_file(os.path.join(doc_dir, "index.html"), file_link_dict, "Historical Spheral Timing Benchmarks", [])
        end_time = time.time()
        run_time = end_time - start_time
        print(f"Script time: {run_time:.4f} seconds")

if __name__=="__main__":
    __spec__ = None
    main()
