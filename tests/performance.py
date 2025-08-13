#!/user/bin/env python3

# This file runs and compares performance tests through the ats system.
# Run using: ./spheral-ats tests/performance.py

import sys, shutil, os, time, stat
import numpy as np
import SpheralConfigs
from SpheralUtilities import TimerMgr
from SpheralTestUtilities import num_3d_cyl_nodes
from ats import configuration

if (not TimerMgr.timers_usable()):
    log("WARNING: Timers not enabled, skipping performance tests. Configure Spheral w/ -DSPHERAL_ENABLE_TIMERS=On", echo=True)
    sys.exit(0)

# Get options from ats
opts = getOptions()

# Adding --threads to the command line arguments of spheral-ats
# can force performance to use multiple threads per rank
num_threads = 1
if "threads" in opts:
    num_threads = opts["threads"]

# Adding --ciRun to the command line arguments of spheral-ats
# triggers copy of Caliper files to benchmark location
benchmark_dir = None
test_runs = 1 # Number of times to run each test
CIRun = False
if "cirun" in opts and opts["cirun"]:
    CIRun = True
    test_runs = 5
    benchmark_dir = opts["benchmark_dir"]
is_rerun = False
if "rerun" in opts and opts["rerun"]:
    is_rerun = True
#---------------------------------------------------------------------------
# Hardware configuration
#---------------------------------------------------------------------------
# This should be {$SYS_TYPE}_{compiler name}_{compiler version}_{mpi or cuda info}
spheral_install_config = SpheralConfigs.config()
mpi_enabled = SpheralConfigs.enable_mpi()
# Retrieve the host name and remove any numbers
temp_uname = os.uname()
hostname = temp_uname[1].rstrip("0123456789")
mac_procs = {"rzhound": 112, "rzwhippet": 112, "ruby": 56,
             "rzadams": 84, "rzvernal": 64, "tioga": 64,
             "rzansel": 40, "lassen": 40, "rzgenie": 36}
# Find out how many nodes our allocation has grabbed
num_nodes = max(1, configuration.machine.numNodes)
if (not mpi_enabled):
    if (num_nodes > 1):
        raise Exception("Should not use more than 1 node when MPI is off")

total_num_cores = 0
try:
    total_num_cores = mac_procs[hostname] * num_nodes
except:
    log("Machine name not recognized", echo=True)
    raise Exception
# If MPI is turned off, thread the whole node
if (not mpi_enabled):
    num_cores = int(total_num_cores)
else:
    # Ideally, tests should be run with 2 nodes and each test will
    # use one entire node, except the 2D tests which use half a node
    num_cores = int(total_num_cores/2)

#---------------------------------------------------------------------------
# Test configurations
#---------------------------------------------------------------------------

# General number of SPH nodes per core
# 5k-10k nodes per core for 3d, 1k nodes per core for 2d
n_per_core_3d = 8000
n_per_core_2d = 1000

Ntotal = int(num_cores*n_per_core_3d)

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
        log(f"Copying {cali_filename} to {outdir}", echo=True)
        if (CIRun):
            shutil.copy(cfile, outfile)
            os.chmod(outfile, perms)
            shutil.chown(outfile, group="sduser")
    if (CIRun):
        cpaths = [outdir, macpath, instpath, benchmark_dir]
        for p in cpaths:
            os.chmod(p, perms)
            shutil.chown(p, group="sduser")

def spheral_setup_test(test_file, test_name, inps, ncores, threads=1, **kwargs):
    '''
    General method for creating an individual performance test
    Parameters:
    test_file: Path to testing script
    test_name: Unique name for test, will be used in Caliper file name
    inps: Command line inputs for test
    ncores: Total number of cores to use for the test, not number of ranks
    threads: Number of threads per rank
    **kwargs: Any additional keyword arguments to pass to ATS tests routine
    '''
    if (not mpi_enabled):
        threads = ncores
        ncores = 1
    for i in range(test_runs):
        if (test_runs > 1):
            cali_name = f"{test_name}_{i}_{int(time.time())}.cali"
        else:
            cali_name = f"{test_name}_{int(time.time())}.cali"
        ccores = int(ncores / threads)
        timer_cmds = f"--caliperFilename {cali_name} --adiakData 'test_name: {test_name}'"
        finps = f"{inps} {timer_cmds}"
        t = test(script=test_file, clas=finps,
                 label=test_name,
                 np=ccores,
                 nt=threads,
                 caliper_filename=cali_name,
                 **kwargs)

# If running CI, run gather_files on exit
if (benchmark_dir):
    onExit(gather_files)
glue(keep=True, independent=True)

#---------------------------------------------------------------------------
# Taylor impact test
#---------------------------------------------------------------------------
test_loc = "functional/Strength/TaylorImpact"

test_file = "TaylorImpact.py"
test_path = os.path.join(test_loc, test_file)
test_name = "3DTAYLOR"

rlen = 0.945
zlen = 7.5
steps = 5
# Estimate nr and nz so the 3D cylindrical node generator creates Ntotal SPH nodes
nz0 = int(np.cbrt(Ntotal)) # Initial guess for nz
nr0 = max(4, int(nz0/4)) # Initial guess for nr
nr, nz = num_3d_cyl_nodes(0., rlen, 0., zlen, 0., 2.*np.pi, nr0, nz0, Ntotal)
gen_inps = f"--geometry 3d --steps {steps} --compatibleEnergy False "+\
    "--clearDirectories False --baseDir None "+\
    "--vizTime None --vizCycle None --siloSnapShotFile None "+\
    f"--rlength {rlen} --zlength {zlen} --nr {nr} --nz {nz}"

# Test variations
test_inp = {"CRK": "--hydroType CRKSPH --densityUpdate SumVoronoiCellDensity",
            "FSI": "--hydroType FSISPH",
            "SOLIDSPH": "--hydroType SPH"}
for tname, tinp in test_inp.items():
    inps = f"{gen_inps} {tinp}"
    spheral_setup_test(test_path, test_name+tname, inps, num_cores, num_threads)

#---------------------------------------------------------------------------
# 3D convection test
#---------------------------------------------------------------------------
test_loc = "unit/Boundary"

test_file = "testPeriodicBoundary-3d.py"
test_path = os.path.join(test_loc, test_file)
test_name = "3DCONV"

# Number of SPH nodes per direction
npd = int(np.cbrt(Ntotal))
steps = 10
inps = f"--nx {npd} --ny {npd} --nz {npd} --steps {steps}"
spheral_setup_test(test_path, test_name, inps, num_cores, num_threads)

#---------------------------------------------------------------------------
# NOH tests
#---------------------------------------------------------------------------
fluid_variations = {"SPH": "--crksph False --solid True",
                    "FSISPH": "--fsisph True --solid True",
                    "CRKSPH": "--crksph True --solid True",
                    "PSPH": "--psph True",
                    "GSPH": "--gsph True",
                    "MFM": "--mfm True",
                    "MFV": "--mfv True"}
test_loc = "functional/Hydro/Noh"

# General input for all Noh tests
gen_noh_inps = "--cfl 0.25 --Cl 1.0 --Cq 1.0 --xfilter 0.0 "+\
    "--nPerh 2.01 --graphics False --clearDirectories False --doCompare False "+\
    "--dataDir None --vizTime None --vizCycle None"

#++++++++++++++++++++
# Noh cylindrical 2d
#++++++++++++++++++++
test_file = "Noh-cylindrical-2d.py"
test_path = os.path.join(test_loc, test_file)
test_name = "NC2D"

steps = 100
rmin = 0.
rmax = 1.
thetaFactor = 0.5
# Only use half the number of cores per node
ncores = int(num_cores/2)
Ntotal2d = n_per_core_2d*ncores
# Determine nRadial to get Ntotal2d number of SPH nodes
# for a constantDTheta distribution
area = np.pi*rmax**2*thetaFactor/2.
dr = np.sqrt(area/Ntotal2d)
nRadial = int(rmax/dr)
gen_inps = f"{gen_noh_inps} --nRadial {nRadial} --steps {steps}"
# Test with different hydro options
for tname, tinp in fluid_variations.items():
    inps = gen_inps + " " + tinp
    spheral_setup_test(test_path, test_name+tname, inps, ncores, num_threads)

#++++++++++++++++++++
# Noh spherical 3d
#++++++++++++++++++++
test_file = "Noh-spherical-3d.py"
test_path = os.path.join(test_loc, test_file)
test_name = "NS3D"

steps = 10
gen_inps = f"{gen_noh_inps} --nx {npd} --ny {npd} --nz {npd} --steps {steps}"
# Test with different hydro options
for tname, tinp in fluid_variations.items():
    inps = gen_inps + " " + tinp
    spheral_setup_test(test_path, test_name+tname, inps, num_cores, num_threads)

# Add a wait to ensure all timer files are done
wait()
