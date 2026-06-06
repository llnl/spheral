#!/user/bin/env python3

# This file runs and compares performance tests through the ats system.
# Run using: ./spheral-ats tests/run_perf.py

import sys, shutil, os
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
test_runs = 1 # Number of times to run each test
if "cirun" in opts and opts["cirun"]:
    test_runs = 3
is_rerun = False
if "rerun" in opts and opts["rerun"]:
    is_rerun = True
#---------------------------------------------------------------------------
# Hardware configuration
#---------------------------------------------------------------------------
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
    '''
    filtered = [test for test in manager.testlist if test.status is PASSED]
    logdir = log.directory
    log(f"Copying caliper files to {logdir}")
    for test in filtered:
        run_dir = test.directory
        cali_filename = test.options["caliper_filename"]
        cfile = os.path.join(run_dir, cali_filename)
        test_name = test.options["label"]
        outfile = os.path.join(logdir, cali_filename)
        shutil.move(cfile, outfile)

# If running CI, run gather_files on exit
onExit(gather_files)
glue(keep=True, independent=True, ngpu=0, nt=1)

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
    inp_dicts = t.create_test_inputs(test_path, test_runs, num_threads)
    for i in inp_dicts:
        test(**i)

# Add a wait to ensure all timer files are done
wait()
