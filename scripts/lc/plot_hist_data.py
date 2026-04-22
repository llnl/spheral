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
# Contains routines for creating the html
import hist_data_utils as hdu

# How many months to plot
num_of_months = 6
# Which region to plot
plot_region = "advance"
subheader = f"""
  <p>Plots show average time per rank per SPH node per step of the {plot_region} region, averaged again over 3 runs.</p>
  <p>3D tests also show the minimum and maximum times of the 3 runs as the bounds of the shaded area.</p>
  <p>Double-click legend entries to isolate them in the plot. Click and drag to zoom.</p>
"""
colors = px.colors.qualitative.G10

def init_worker():
    global reg_names, spot_metric
    reg_names = ["advance", "ConnectivityMap_computeConnectivity", "computeVoronoiVolume"]
    spot_metric = "avg#inclusive#sum#time.duration"

def normalize_time(gls):
    "Normalize time by the number of nodes and number of steps"
    sph_nodes = int(gls['total_internal_nodes'])
    steps = int(gls['total_steps'])
    return 1./(sph_nodes*steps)

def extract_data(cali_file):
    global reg_names, spot_metric
    (records, gls) = cr.read_caliper_contents(cali_file)
    runday = int(gls["launchdate"])
    time_norm = normalize_time(gls)
    rtimes = {}
    for rec in records:
        if ('path' in rec):
            for rg in reg_names:
                if (rec["path"][-1] == rg):
                    tval = float(rec[spot_metric])*time_norm
                    if (rg in rtimes.keys()):
                        rtimes[rg] += tval
                    else:
                        rtimes[rg] = tval
    return runday, rtimes

def get_spec(cali_file):
    # If cali_file is just a dir, grab the latest caliper file
    if os.path.isdir(cali_file):
        perf_dir = sorted(glob.glob(os.path.join(cali_file, "*")))[-1]
        cali_file = glob.glob(os.path.join(perf_dir, "*.cali"))[0]
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

def split_test_name(test_name, test_base_names):
    "Split the test name into the test base and the test variant"
    for base in test_base_names:
        if test_name.startswith(base):
            variant = test_name[len(base):]
            return base, variant
    return None, None

def update_spec(spec):
    "Spec names changes from spack 0.12 to 1.0, so we must update the names to match"
    newspec = spec.replace("clang", "llvm")
    newspec = newspec.replace("rocmcc", "llvm-amdgpu")
    return newspec

# Convert array of dictionaries into panda dataframe
def convert_to_dataframe(in_array, test_base_names):
    data = []
    for in_data in in_array:
        mac_name = in_data["machine"]
        spec_name = update_spec(in_data["spec"])
        config = in_data["config"]
        data_dict = in_data["data_dict"]
        for test_name, tdata in data_dict.items():
            test_base, test_var = split_test_name(test_name, test_base_names)
            dates = tdata["dates"]
            times = tdata["times"]
            for date, time_dict in zip(dates, times):
                for reg, val in time_dict.items():
                    data.append({"Machine": mac_name,
                                 "Spec": spec_name,
                                 "Config": config,
                                 "Test Name": test_name,
                                 "Test Base": test_base,
                                 "Test Var": test_var,
                                 "Region": reg,
                                 "Date": str(date),
                                 "Time": val})
    df = pd.DataFrame(data)
    df = df.sort_values(by="Date")
    return df

def create_plot(in_df):
    mac_name = in_df["Machine"].unique()[0]
    specs = in_df["Spec"].unique()
    test_bases = sorted(in_df["Test Base"].unique())
    figs = []
    # Pick the region to plot
    df = in_df[in_df["Region"] == plot_region]
    for spec in specs:
        spec_df = df[df["Spec"] == spec]
        for test_base in test_bases:
            test_df = spec_df[spec_df["Test Base"] == test_base]
            title = f"Test: {test_base}, Spec: {spec}, Avg"
            if "2D" not in test_base:
                title += ", Min, and Max"
            fig = px.line(title=title)
            test_vars = sorted(test_df["Test Var"].unique())
            for cindx, var in enumerate(test_vars):
                ccolor = colors[cindx]
                cdf = test_df[test_df["Test Var"] == var]
                # Plot the mean for any runs on a given day
                mean_df = cdf.groupby("Date")["Time"].mean().reset_index()
                # Plot the min and max for non-2D runs
                if ("2D" not in test_base):
                    min_df = cdf.groupby("Date")["Time"].min().reset_index()
                    max_df = cdf.groupby("Date")["Time"].max().reset_index()
                    fig.add_scatter(x=max_df["Date"],
                                    y=max_df["Time"],
                                    legendgroup=var,
                                    showlegend=False,
                                    line=dict(color=ccolor))
                    fig.add_scatter(x=min_df["Date"],
                                    y=min_df["Time"],
                                    legendgroup=var,
                                    showlegend=False,
                                    fill="tonexty",
                                    line=dict(color=ccolor))
                fig.add_scatter(x=mean_df["Date"],
                                y=mean_df["Time"],
                                line=dict(color=ccolor),
                                marker=dict(color=ccolor, size=10),
                                mode="lines+markers",
                                name=f"{var}",
                                legendgroup=var,
                                showlegend=True)
                fig.update_yaxes(title_text=f"Grind time of region {plot_region}")
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
    if (rank == 0):
        if (os.path.exists(doc_dir)):
            shutil.rmtree(doc_dir)
        os.makedirs(doc_dir)
    perf_paths = ["tests", "spheral/tests"]
    sys.path.extend([os.path.join(cur_dir, i) for i in perf_paths])
    import perf_tests as pt
    if (args.pkg_name == "LLNLSpheral"):
        import llnl_perf_tests
    test_names = pt.get_all_test_names()
    if (args.test):
        test_names = [test_names[0]]
    test_bases = pt.get_test_bases()
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
            if (len(test_files) > 0):
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
        df_data = convert_to_dataframe(all_data, test_bases)
        macs = df_data["Machine"].unique()
        file_link_dict = {mac: f"{mac}.html" for mac in macs}
        file_link_dict.update({"Home": "index.html"})
        for mac in macs:
            mac_df = df_data[df_data["Machine"] == mac]
            figs = create_plot(mac_df)
            page_content = hdu.create_page_content(figs)
            hdu.create_html_file(file_name=os.path.join(doc_dir, f"{mac}.html"),
                                 file_dict=file_link_dict,
                                 title=mac,
                                 content=page_content,
                                 subheader=subheader)
        hdu.create_html_file(file_name=os.path.join(doc_dir, "index.html"),
                             file_dict=file_link_dict,
                             title="Historical Spheral Timing Benchmarks")
        end_time = time.time()
        run_time = end_time - start_time
        print(f"Script time: {run_time:.4f} seconds")

if __name__=="__main__":
    __spec__ = None
    main()
