#-------------------------------------------------------------------------------
# Functions for adding Caliper and Adiak parsing arguments and initializing
# the timer manager
#-------------------------------------------------------------------------------

import argparse, mpi
import SpheralConfigs
from SpheralUtilities import TimerMgr
from SpheralUtilities import adiak_value
import SpheralOpenMP

cali_args = ["caliperConfig", "caliperFilename", "caliperConfigJSON"]

def parse_dict(string):
    """
    Function to parse a dictionary provided through the command line
    """
    try:
        inp_dict = dict(item.split(":") for item in string.split(","))
    except:
        raise SyntaxError("Input to --adiakData must be in key:value format, separated by commas")
    new_dict = {}
    for ikey, ival in inp_dict.items():
        try:
            key = eval(ikey)
        except:
            key = ikey.strip()
        try:
            val = eval(ival)
        except:
            val = ival.strip()
        new_dict.update({key: val})
    return new_dict

def add_timing_args(parser):
    """
    Add Caliper and Adiak arguments to the parser
    """
    # Allow Adiak values to be set on the command line
    # Inputs are a string that can be evaluated into a dictionary
    # For example, --adiakData "testname: ShockTube1, testing:3"
    parser.add_argument("--adiakData", default=None, action='append',
                        type=parse_dict)
    # This logic checks if the user already set a Caliper
    # argument and default value and prevents adding the argument
    # if it already exists
    arg_list = [action.dest for action in parser._actions]
    for ca in cali_args:
        if (ca not in arg_list):
            parser.add_argument(f"--{ca}", default="", type=str)

def init_timer(args):
    """
    Initializes the timing manager and adds input values to Adiak from parsed arguments
    """
    if args.verbose:
        if (args.caliperConfig):
            print("  *  caliperConfig = ", args.caliperConfig)
        if (args.caliperFilename):
            print("  *  caliperFilename = ", ars.caliperFilename)
        if (args.caliperConfigJSON):
            print("  *  caliperConfigJSON = ", args.caliperConfigJSON)
    if (not TimerMgr.timers_usable()):
        if (args.caliperConfig or args.caliperFilename or args.caliperConfigJSON):
            print("WARNING: Caliper command line inputs provided for "+\
                  "non-timer install. Reconfigure the install with "+\
                  "-DSPHERAL_ENABLE_TIMERS=ON to be able to use Caliper timers.")
    if(args.caliperConfigJSON):
        TimerMgr.load(args.caliperConfigJSON)
        if(not args.caliperConfig):
            raise RuntimeError("SpheralOptionParser: specifying a configuration file without "+\
                               "using one of the configurations means no timers are started")
    off_tests = ["none", "off", "disable", "disabled", "0"]
    # Check if Caliper is turned off
    if (args.caliperConfig):
        if (args.caliperConfig.lower() in off_tests):
            return
        TimerMgr.add(args.caliperConfig)
        TimerMgr.start()
    else:
        import os, sys
        # If output name for Caliper is given, use it
        if (args.caliperFilename):
            testname = args.caliperFilename
        else:
            from datetime import datetime
            # Append the current day and time to the filename
            unique_digits = datetime.now().strftime("_%Y_%m_%d_%H%M%S_%f")
            # Name file based on name of python file being run
            testname = os.path.splitext(os.path.basename(sys.argv[0]))[0]
            testname += unique_digits + ".cali"
        TimerMgr.default_start(testname)
    # Add number of ranks and threads per rank
    adiak_value("threads_per_rank", SpheralOpenMP.omp_get_num_threads())
    adiak_value("num_ranks", mpi.procs)
    # Add metadata about install
    adiak_value("install_config", SpheralConfigs.config())
    adiak_value("build_type", SpheralConfigs.build_type())
    adiak_value("spec", SpheralConfigs.spec())
    if SpheralConfigs.git_hash():
        adiak_value("git_hash", SpheralConfigs.git_hash())
        adiak_value("git_branch", SpheralConfigs.git_branch())

    # Add --adiakData command line inputs as Adiak metadata
    if (args.adiakData):
        for i in args.adiakData:
            for key, val in i.items():
                adiak_value(key, val)
    return
