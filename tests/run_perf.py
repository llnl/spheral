#!/user/bin/env python3

# This file runs and compares performance tests through the ats system.
# Run using: ./spheral-ats tests/run_perf.py

import sys, shutil, os, time, stat
import numpy as np
import SpheralConfigs
from SpheralTestUtilities import num_3d_cyl_nodes
from ats import configuration

if (not SpheralConfigs.timers_enabled()):
    log("WARNING: Timers not enabled, skipping performance tests. Configure Spheral w/ -DSPHERAL_ENABLE_TIMERS=On", echo=True)
    sys.exit(0)

# Get options from ats
opts = getOptions()

# Get install test path
test_path = os.path.join(opts["install_path"], "tests")

# Adding --threads to the command line arguments of spheral-ats
# can force performance to use multiple threads per rank
num_threads = 1
if "threads" in opts:
    num_threads = opts["threads"]

# Adding --ciRun to the command line arguments of spheral-ats
# triggers move of Caliper files to benchmark location
benchmark_dir = None
test_runs = 1 # Number of times to run each test
CIRun = False
if "cirun" in opts and opts["cirun"]:
    CIRun = True
    test_runs = 3
    benchmark_dir = opts["benchmark_dir"]
is_rerun = False
if "rerun" in opts and opts["rerun"]:
    is_rerun = True
#---------------------------------------------------------------------------
# Hardware configuration
#---------------------------------------------------------------------------
# This should be {$SYS_TYPE}_{compiler name}_{compiler version}_{mpi or cuda info}
spheral_install_config = SpheralConfigs.config()
mpi_enabled = SpheralConfigs.mpi_enabled()
# Retrieve the host name and remove any numbers
temp_uname = os.uname()
hostname = temp_uname[1].rstrip("0123456789")
mac_procs = {"rzhound": 112, "rzwhippet": 112, "dane": 112,
             "rzadams": 96, "rzvernal": 64, "tioga": 64,
             "rzgenie": 36}
# Find out how many nodes our allocation has grabbed
num_nodes = max(1, configuration.machine.numNodes)
if (not mpi_enabled):
    if (num_nodes > 1):
        raise Exception("Should not use more than 1 node when MPI is off")

num_cores = 0
try:
    num_cores = mac_procs[hostname] * num_nodes
except:
    log("Machine name not recognized", echo=True)
    raise Exception

def gather_files(manager):
    '''
    Function to gather Caliper file when ATS is finished running.
    Used by ATS for gathering benchmark Caliper files.
    '''
    instpath = os.path.join(benchmark_dir, spheral_install_config)
    macpath = os.path.join(instpath, hostname)
    outdir = os.path.join(macpath, "latest")
    if (not is_rerun):
        if (os.path.exists(outdir)):
            # Move existing benchmark data to a different directory
            log(f"Renaming existing {outdir} to {int(time.time())}", echo=True)
            os.rename(outdir, os.path.join(macpath, f"{int(time.time())}"))
        log(f"Creating {outdir}", echo=True)
        os.makedirs(outdir)
    filtered = [test for test in manager.testlist if test.status is PASSED]
    # Set read/write/execute permissions for owner and group
    perms = stat.S_IRWXU | stat.S_IRWXG
    for test in filtered:
        run_dir = test.directory
        cali_filename = test.options["caliper_filename"]
        cfile = os.path.join(run_dir, cali_filename)
        test_name = test.options["label"]
        outfile = os.path.join(outdir, cali_filename)
        log(f"Moving {cali_filename} to {outdir}", echo=True)
        if (CIRun):
            shutil.move(cfile, outfile)
            os.chmod(outfile, perms)
            shutil.chown(outfile, group="sduser")
    if (CIRun):
        cpaths = [outdir, macpath, instpath, benchmark_dir]
        for p in cpaths:
            os.chmod(p, perms)
            shutil.chown(p, group="sduser")

def spheral_setup_test(test_param, threads=1, **kwargs):
    '''
    General method for creating an individual performance test
    Parameters:
    test_param: TestParam class for specific test
    ncores: Total number of cores to use for the test, not number of ranks
    threads: Number of threads per rank
    **kwargs: Any additional keyword arguments to pass to ATS tests routine
    '''
    ncores = test_param.ncores
    if (not mpi_enabled):
        threads = ncores
        ncores = 1
    test_file = test_param.get_test_file(test_path)
    tests = test_param.get_tests()
    for test_name, inps in tests.items():
        for i in range(test_runs):
            if (test_runs > 1):
                cali_name = f"{test_name}_{i}_{int(time.time())}.cali"
            else:
                cali_name = f"{test_name}_{int(time.time())}.cali"
            timer_cmds = f"--caliperFilename {cali_name} --adiakData 'test_name: {test_name}'"
            finps = f"{inps} {timer_cmds}"
            t = test(script=test_file, clas=finps,
                     label=test_name,
                     np=ncores,
                     nt=threads,
                     caliper_filename=cali_name,
                     **kwargs)

# If running CI, run gather_files on exit
if (benchmark_dir):
    onExit(gather_files)
glue(keep=True, independent=True)

#---------------------------------------------------------------------------
# Test configurations
#---------------------------------------------------------------------------
perf_paths = ["", "spheral/tests", "llnlspheral/tests"]
sys.path.extend([os.path.join(test_path, i) for i in perf_paths])
import perf_tests as pt
if(os.path.exists(perf_paths[-1])):
    import llnl_perf_tests

all_tests = pt.get_all_tests(num_cores)

for t in all_tests:
    spheral_setup_test(t, num_threads)

# Add a wait to ensure all timer files are done
wait()
