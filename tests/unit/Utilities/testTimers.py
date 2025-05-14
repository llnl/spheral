#
#
#ATS:test(SELF, "--caliperFilename 'timer_test_1.cali'", label="Timer test 1", np=1)
#ATS:test(SELF, "--caliperConfig 'none'", label="Timer test 2", np=1)
#ATS:test(SELF, "--caliperFilename 'timer_test_3.cali' --adiakData 'adiak_test: 1, test_adiak: two'", label="Timer test 3", np=1)
#ATS:test(SELF, "--caliperFilename 'timer_test_4.cali' --adiakData 'adiak_test: 1, test_adiak: two'", label="Timer test 3", np=8)
#

import Spheral
from SpheralTestUtilities import *
from SpheralOptionParser import *
from SpheralUtilities import *
import mpi
import SpheralConfigs

import sys, os, time

# Test set Adiak inputs
test_dict = {"perf_test": "weak_scaling",
             "rank_count": mpi.procs,
             "fake_float": 2.141}
for key, val in test_dict.items():
    adiak_value(key, val)

# Test the --adiakData input. This must match what is
# hard-coded in the ATS magic lines
adiak_data_dict = {"adiak_test": 1, "test_adiak": "two"}

commandLine()

# Remove cali files from previous test runs
caliper_file = TimerMgr.get_filename()
if (os.path.exists(caliper_file)):
    if (mpi.rank == 0):
        os.remove(caliper_file)

do_timers = False
if (TimerMgr.is_started()):
    do_timers = True

run_count = 8
sleep_time = 1.E-4
fake_timer_name = "test_timer"

for i in range(run_count):
    TimerMgr.timer_begin(fake_timer_name)
    time.sleep(sleep_time)
    TimerMgr.timer_end(fake_timer_name)

# Read in Caliper file and process it
if (do_timers and TimerMgr.get_filename()):
    adiak_fini()
    TimerMgr.fini()
    mpi.barrier()
    caliper_loc = SpheralConfigs.caliper_module_path()
    if (not caliper_loc):
        raise FileNotFoundError("Caliper file not found")
    sys.path.append(caliper_loc)
    import caliperreader as cr
    if (mpi.rank == 0):
        r = cr.CaliperReader()
        r.read(caliper_file)
        records = r.records

        # Test for timer name
        assert fake_timer_name in records[1]['region'], f"{fake_timer_name} timer not found"

        # Test for function count
        count_val = int(eval(records[1]["avg#sum#rc.count"]))
        assert count_val == run_count, "Caliper function count is off"

        # Note: CaliperReader reads everything as strings for some terrible reason
        # we must convert the Adiak values first
        adiak_inp = {}
        for key, val in r.globals.items():
            try:
                newval = eval(val)
            except:
                newval = val
            adiak_inp.update({key: newval})

        # Test Adiak output for explicitly set values
        assert test_dict.items() <= adiak_inp.items(),\
            "incorrect Adiak values found in Caliper file"

        # Test --adiakData command line input
        if ("adiakData" in adiak_inp):
            assert adiak_data_dict.items() <= adiak_inp.items(),\
                "incorrect adiakData inputs found in Caliper file Adiak values"
