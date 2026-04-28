"""
Compare performance data for Spheral

Do a performance regression test on LC systems with the following steps:

  1. Run the performance test in Spheral

     $> ./spheral-ats --numNodes 2 --logs test_dir_name tests/performance.py

  2. Run this script and point to the directory created by ATS in step 1

     $> ./spheral performance_analysis.py --perfdata test_dir_name
"""

import os, sys, shutil, glob
import argparse

try:
    import thicket as th
    import hatchet as ht
except:
    print("Thicket not found. Make sure the TPLs are up-to-date.")
    raise Exception

from IPython.display import display
from IPython.display import HTML

import numpy as np

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

# Set the region and metric to use for comparisons
comp_region = "advance"
comp_metric = "Avg time/rank" # Inclusive timer
# Set the metric to use for display
disp_metric = "Avg time/rank (exc)" # Exclusive timer

percent = 0.08
def compute_threshold(sf, metric):
    """
    Compute the threshold of time as thresh = 0.08*mu + 2*sigma
    Will be used in the comparison abs(t - mu) > thresh to determine
    notable differences in times
    Parameters:
    sf: Thicket statsframe (Hatchet Graphframe) for reference data
    metric: Metric for threshold
    Returns:
    thresh: Time limit threshold Graphframe variable
    """
    mu = sf.dataframe[metric+"_mean"]
    if (metric+"_std" in sf.dataframe):
        sigma = sf.dataframe[metric+"_std"]
        return percent*mu + 2.*sigma
    return percent*mu

def check_for_region(data, region):
    "Check if region exists in Thicket/Hatchet"
    return not data.dataframe[data.dataframe["name"] == region].empty

def get_times(gf, region, metric = "Avg time/rank"):
    """
    Return the times for a given region and metric.
    If multiple regions have the same name, the sum of those times are returned.
    Returns an array of times for each profile in the Thicket
    """
    if (type(gf) is th.thicket.Thicket):
        cregion = gf.get_node(region)
        profs = gf.profile
        if len(profs) > 1:
            times = []
            for i in profs:
                vth = gf.filter_profile([i])
                times.append(sum(vth.dataframe.loc[cregion][metric].values))
        else:
            times = [sum(gf.dataframe.loc[cregion][metric].values)]
    else:
        names = gf.dataframe["name"]
        indx = names[names == region].index
        times = [sum(gf.dataframe.loc[indx][metric].values)]
    return times

def remove_nans(gf, metric="Avg time/rank"):
    "Remove rows with NANs in a GraphFrame or Thicket or list/dict of those types"
    if (type(gf) is dict or type(gf) is th.groupby.GroupBy):
        newdict = {}
        for key, val in gf.items():
            newval = remove_nans(val)
            newdict.update({key: newval})
        return newdict
    elif (type(gf) is list):
        newlist = []
        for val in gf:
            newval = remove_nans(val)
            newlist.append(newval)
        return newlist
    elif (type(gf) is th.thicket.Thicket):
        query = th.query.Query().match("+", lambda row: row[metric].apply(lambda x: not np.isnan(x)).all())
        return gf.query(query)
    elif (type(gf) is ht.GraphFrame):
        query = ht.query.Query().match("+", lambda x: not np.isnan(x[metric]))
        return gf.filter(query)
    else:
        raise TypeError(f"Unrecognized type in remove_nans {type(gf)}")

def group_tests(data):
    """
    Groups input data based on tests.
    Parameters:
    data: Thicket to filter
    Returns:
    filt: GroupBy of Thickets based on tests
    """
    test_group = ["test_name", "total_internal_nodes", "total_steps"]
    return data.groupby(test_group)

def compare_metadata(cdata, rdata, tests):
    cmdata = cdata.get_unique_metadata()
    rmdata = rdata.get_unique_metadata()
    failed_configs = []
    for t in tests:
        cval = cmdata[t]
        rval = rmdata[t]
        if (cval != rval):
            failed_configs.append(f"{t}: {cval} != {rval}")
    return failed_configs

def compare_tests(cdata, rdata):
    "Determine why a test configuration differs"
    tests = ["total_internal_nodes", "total_steps"]
    return compare_metadata(cdata, rdata, tests)

def compare_config(cdata, rdata):
    """
    Check if tests between two sets of data are compatible for direct comparison.
    Do not use if comparing across hardware or install configurations.
    Assumes Thicket data has comparable test_name, total_internal_nodes, and total_steps
    """
    hardware_tests = ["install_config", "cluster", "jobsize", "threads_per_rank"]
    return compare_metadata(cdata, rdata, hardware_tests)

def filter_tests(data, test_name):
    return data.filter_metadata(lambda x: x["test_name"] == test_name)

def get_caliper_files_and_bench(file_path):
    atsFile = os.path.join(file_path, "atsr.py")
    cali_files = []
    benchmarks = None
    # If caliper file is provided directly
    if (".cali" in file_path):
        cali_files = [file_path]
    # If perf-dir is an ATS output directory, find the Caliper files from atsr.py
    elif (os.path.exists(atsFile)):
        # Run atsr.py and put values into globals
        exec(compile(open(atsFile).read(), atsFile, 'exec'), globals())
        state = globals()["state"]
        tests = [t for t in state["testlist"] if t['status'] == PASSED]
        for test in tests:
            if ("caliper_filename" in test["options"]):
                cali_file = test["options"]["caliper_filename"]
            else:
                raise RuntimeError("This tool only works on ATS output from run_perf.py")
            # Check if benchmark_dir is in ats options
            if (not benchmarks and "benchmark_dir" in test["options"]):
                benchmarks = test["options"]["benchmark_dir"]
            # Check if the run_perf.py gathered the Caliper files or not
            cfile = os.path.join(file_path, cali_file)
            if (not os.path.exists(cfile)):
                cfile = os.path.join(test["directory"], cali_file)
            cali_files.append(cfile)
    else:
        newpath = os.path.join(file_path, "**/*.cali")
        print(f"Searching {newpath}")
        cali_files = glob.glob(newpath, recursive=True)
        if (not cali_files):
            raise RuntimeError(f"No caliper files found in {file_path}")
    return cali_files, benchmarks

def get_caliper_files(file_path):
    cali_files, unused_bench = get_caliper_files_and_bench(file_path)
    return cali_files

def main():
    #---------------------------------------------------------------------------
    # Setup argument parser
    #---------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        usage="""If doing a performance regression test, use inputs:
        --perfdata1 /path/to/perfdata --ref /path/to/benchmark/data
        Otherwise, compare two performance outputs using:
        --perfdata1 /path/to/perfdata --perfdata2 /path/to/perfdata2
        If only --perfdata1 is specified, --ref is set to be the latest upstream benchmark data
        If only --perfdata is specified, a Thicket tree is displayed of the timers.
        """)
    group1 = parser.add_mutually_exclusive_group()
    group1.add_argument("--perfdata1", type=str, default=None,
                        help="Directory containing an atsr.py file "+\
                        "or a collection of Caliper files. Use when doing comparisons "+\
                        "with --perfdata2 or --ref/benchmark data")
    group1.add_argument("--perfdata", type=str, default=None,
                        help="Directory containing an atsr.py file "+\
                        "or a collection of Caliper files.")
    group2 = parser.add_mutually_exclusive_group()
    group2.add_argument("--ref", type=str, default=None,
                       help="Directory of Caliper files to use as reference for "+\
                       " comparing against a threshold. Exits with failure perfdata1 "+\
                       "times exceed threshold.")
    group2.add_argument("--perfdata2", type=str, default=None,
                       help="Directory of an atsr.py file or a collection of Caliper "+\
                       "files to compare against perfdata1.")
    parser.add_argument("--test-name", type=str, default=None,
                        help="If comparing a specific test, default is compare all.")
    parser.add_argument("--display", action="store_true",
                        help="Display a tree for timers that failed.")
    args = parser.parse_args()

    # Create a Thicket of the current performance data
    #-------------------------------------------------
    if (args.perfdata1):
        perfdata = args.perfdata1
        no_comp = False
    elif (args.perfdata):
        perfdata = args.perfdata
        no_comp = True
    else:
        raise Exception("Must specify either --perfdata or --perfdata1")
    if (not os.path.exists(perfdata)):
        raise Exception(f"Cannot find {perfdata}")
    cali_files, benchmarks = get_caliper_files_and_bench(perfdata)
    if (len(cali_files) == 0):
        raise Exception(f"No .cali files found in {perfdata}")
    curdata = th.Thicket.from_caliperreader(cali_files, disable_tqdm=True)
    # Filter data set by tests
    cur_test_data = group_tests(curdata)
    cur_test_data = remove_nans(cur_test_data)
    if (no_comp):
        if (args.perfdata2 or args.ref):
            raise RuntimeError("To compare with --perfdata2 or --ref, use --perfdata1")
        for test_key, ctest in cur_test_data.items():
            metric = "Avg time/rank"
            if (len(ctest.profile) > 1):
                th.stats.mean(ctest, [metric])
                metric += "_mean"
            else:
                ctest.move_metrics_to_statsframe([metric])
            display(ctest.statsframe.tree(metric))
        sys.exit()

    # Create a Thicket of the other performance data
    #-----------------------------------------------
    do_thresh_test = False
    if (args.perfdata2):
        if (not os.path.exists(args.perfdata2)):
            raise Exception(f"--perfdata2 location {args.perfdata2} does not exist")
        cali_ref_files = get_caliper_files(args.perfdata2)
    else:
        do_thresh_test = True
        if (args.ref):
            if (not os.path.exists(args.ref)):
                raise Exception(f"--ref location {args.ref} does not exist")
            cali_ref_files = get_caliper_files(args.ref)
        else:
            # If no ref or benchmark_dir is provided, look for the benchmark in the
            # atsr.py file
            ref_files = benchmarks
            # Check in a Caliper file
            if (not ref_files):
                try:
                    ref_files = curdata.metadata["benchmark_dir"].iloc[0]
                except:
                    raise Exception("No reference or benchmark data specified")
            # If using benchmark reference data, only grab the current install/machine
            # Get install config and machine name from current data
            install_config = curdata.metadata["install_config"].iloc[0]
            machine_name = curdata.metadata["cluster"].iloc[0]
            ref_loc = os.path.join(ref_files, install_config, machine_name)
            if (not os.path.exists(ref_loc)):
                raise Exception(f"Benchmark location {ref_loc} does not exists")
            cali_ref_files = glob.glob(os.path.join(ref_loc, "**/*.cali"), recursive=True)

    if (len(cali_ref_files) == 0):
        raise Exception(f"No Caliper files found in {cali_ref_files}")

    test_status = {}
    # Iterate over each test
    for test_key, ctest in cur_test_data.items():
        test_name = test_key[0]
        if (args.test_name and args.test_name != test_name):
            continue
        ref_test_files = [i for i in cali_ref_files if test_name in i]
        ref_test_data = th.Thicket.from_caliperreader(ref_test_files, disable_tqdm=True)
        ref_test_data = group_tests(ref_test_data)
        test_sph_nodes = test_key[1]
        test_steps = test_key[2]
        if (test_key not in ref_test_data):
            # This means the test configurations differ (number of time steps etc)
            if (test_name not in refdata.get_unique_metadata()["test_name"]):
                skip_msg = f"{test_name} not found in reference data"
                test_status.update({test_name: ("SKIPPED-TEST", [skip_msg])})
            else:
                rtest = filter_tests(refdata, test_name)
                ftest_configs = compare_tests(ctest, rtest)
                test_status.update({test_name: ("SKIPPED-TEST", ftest_configs)})
            continue
        rtest = ref_test_data[test_key]
        if (do_thresh_test):
            fh_configs = compare_config(ctest, rtest)
            if (fh_configs):
                # This means the hardware/compiler configurations differs
                test_status.update({test_name: ("SKIPPED-CONF", fh_configs)})
                continue
        metrics = [comp_metric, disp_metric]
        cmetrics = [x+"_mean" for x in metrics]
        # Get stats for current tests
        th.stats.mean(ctest, metrics)
        th.stats.mean(rtest, metrics)
        # Get stats for other tests
        if (len(rtest.profile) > 1):
            th.stats.std(rtest, metrics)
        # Extract times of comp_region
        if (not check_for_region(ctest, comp_region)):
            print(f"{comp_region} not found in {perfdata}")
            continue
        if (not check_for_region(rtest, comp_region)):
            print(f"{comp_region} not found in {os.path.dirname(ref_test_files[0])}")
            continue
        cmain = get_times(ctest.statsframe, comp_region, cmetrics[0])[0]
        rmain = get_times(rtest.statsframe, comp_region, cmetrics[0])[0]
        main_diff = cmain - rmain
        if (do_thresh_test):
            # Compute the max allowable time for the comp_region
            ctest.statsframe.dataframe["thresh"] = compute_threshold(rtest.statsframe, metrics[0])
            ref_thresh = get_times(ctest.statsframe, comp_region, "thresh")[0]
            if (main_diff > ref_thresh):
                cur_status = "FAILED"
                if args.display:
                    # Display the relative difference of the exclusive avg time/rank
                    vals1 = ctest.statsframe.dataframe[cmetrics[1]]
                    vals2 = rtest.statsframe.dataframe[cmetrics[1]]
                    ctest.statsframe.dataframe["exc_rel_diff_percent"] = (vals1/vals2 - 1.)*100.
                    display(ctest.statsframe.tree("exc_rel_diff_percent", cmetrics[1]))
            elif (main_diff < -ref_thresh):
                cur_status = "PASSED"
                if args.display:
                    vals2 = rtest.statsframe.dataframe[cmetrics[1]]
                    ctest.statsframe.dataframe["pdata2"] = vals2
                    display(ctest.statsframe.tree(cmetrics[1], "pdata2"))
            else:
                cur_status = "PASSED"
            test_status.update({test_name: (cur_status, cmain, rmain, ref_thresh)})
        else:
            test_status.update({test_name: ("NA", cmain, rmain)})
            if args.display:
                ctest.statsframe.dataframe["pdata2"] = rtest.statsframe.dataframe[cmetrics[1]]
                display(ctest.statsframe.tree(cmetrics[1], "pdata2"))
    if (do_thresh_test):
        print("Negative values mean local data was faster than reference")
        print(f"Test name: test status, % change in time of {comp_region} region")
    else:
        print("Negative values mean perfdata1 was faster than perfdata2")
        print(f"Test name: % change in time of {comp_region} region")
    failed_tests = []
    for test_name, val in test_status.items():
        if ("SKIPPED" in val[0]):
            diff_str = " ".join(str(x) for x in val[1])
            print(f"{test_name}: SKIPPED, Differences found in: {diff_str}")
        elif (do_thresh_test):
            ctime = val[1]
            rtime = val[2]
            thresh = val[3]
            if ("FAILED" in val[0]):
                failed_tests.append(test_name)
                print(f"{test_name}: FAILED, {(ctime/rtime-1.)*100.:0.3f}%")
            else:
                print(f"{test_name}: PASSED, {(ctime/rtime-1.)*100.:0.3f}%")
        else:
            ctime = val[1]
            rtime = val[2]
            print(f"{test_name}: {(ctime/rtime-1.)*100.:0.3f}%")
    if (len(failed_tests) > 0):
        print("ERROR: The following tests have failed:")
        for i in failed_tests:
            print(i)

if __name__ == "__main__":
    main()
