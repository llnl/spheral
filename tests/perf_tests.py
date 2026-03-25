#!/user/bin/env python3

# Configure the tests to be run by run_perf.py

import os
import numpy as np

# General number of SPH nodes per core
# 5k-10k nodes per core for 3d, 1k nodes per core for 2d
n_per_core_3d = 8000
n_per_core_2d = 1000

# Template class for tests. Any class that inherits from this represents a set of tests
class TestParams:
    # Constructor
    def __init__(self, test_name, test_file, test_vars = None):
        self.test_type = "spheral"
        self.test_name = test_name
        self.test_file = test_file
        self.test_vars = test_vars
        self.ncores = None
        self.gen_inp = None

    def get_test_file(self, install_test_path):
        # Determine path to test files either tests/ or tests/spheral/tests
        test_path_1 = os.path.join(install_test_path, self.test_type, "tests", self.test_file)
        test_path_2 = os.path.join(install_test_path, self.test_file)
        if(os.path.exists(test_path_1)):
            return test_path_1
        elif(os.path.exists(test_path_2)):
            return test_path_2
        else:
            raise OSError(f"Test file {self.test_file} cannot be found.")

    def test_names(self):
        if not self.test_vars:
            return [self.test_name]
        else:
            return [self.test_name + i for i in list(self.test_vars.keys())]

    def set_gen_inputs(self):
        pass

    def get_tests(self):
        self.set_gen_inputs()
        if not self.test_vars:
            return {self.test_name: self.gen_inp}
        else:
            return {f"{self.test_name}{tvar}": f"{tval} {self.gen_inp}" for tvar, tval in self.test_vars.items()}

#---------------------------------------------------------------------------
# Taylor impact test
#---------------------------------------------------------------------------        
class TaylorImpact(TestParams):
    def __init__(self, ncores):
        super().__init__("3DTAYLOR",
                         "functional/Strength/TaylorImpact/TaylorImpact.py",
                         {"CRK": "--hydroType CRKSPH --densityUpdate SumVoronoiCellDensity",
                          "FSI": "--hydroType FSISPH",
                          "SOLIDSPH": "--hydroType SPH"})
        # Only use half the number of cores
        self.ncores = int(ncores/2)

    def set_gen_inputs(self):
        from SpheralTestUtilities import num_3d_cyl_nodes
        Ntotal = self.ncores*n_per_core_3d
        rlen = 0.945
        zlen = 7.5
        steps = 5
        # Estimate nr and nz so the 3D cylindrical node generator creates Ntotal SPH nodes
        nz0 = int(np.cbrt(Ntotal)) # Initial guess for nz
        nr0 = max(4, int(nz0/4)) # Initial guess for nr
        nr, nz = num_3d_cyl_nodes(0., rlen, 0., zlen, 0., 2.*np.pi, nr0, nz0, Ntotal)
        self.gen_inp = f"--geometry 3d --steps {steps} --compatibleEnergy False "+\
            "--clearDirectories False --baseDir None "+\
            "--vizTime None --vizCycle None --siloSnapShotFile None "+\
            f"--rlength {rlen} --zlength {zlen} --nr {nr} --nz {nz}"

#---------------------------------------------------------------------------
# 3D convection test
#---------------------------------------------------------------------------
class Conv3D(TestParams):
    def __init__(self, ncores):
        super().__init__("3DCONV",
                         "unit/Boundary/testPeriodicBoundary-3d.py")
        # Only use half the total cores
        self.ncores = int(ncores/2)

    def set_gen_inputs(self):
        Ntotal = self.ncores*n_per_core_3d
        npd = int(np.cbrt(Ntotal))
        steps = 10
        self.gen_inp = f"--nx {npd} --ny {npd} --nz {npd} --steps {steps}"

#---------------------------------------------------------------------------
# 2D NOH tests
#---------------------------------------------------------------------------
class NOH2D(TestParams):
    def __init__(self, ncores):
        super().__init__("NC2D",
                         "functional/Hydro/Noh/Noh-cylindrical-2d.py",
                         {"SPH": "--crksph False --solid True",
                          "FSISPH": "--fsisph True --solid True",
                          "CRKSPH": "--crksph True --solid True",
                          "PSPH": "--psph True",
                          "GSPH": "--gsph True",
                          "MFM": "--mfm True",
                          "MFV": "--mfv True"})
        # Only use 1/4 the number of cores
        self.ncores = int(ncores/4)
        self.noh_gen_inp = "--cfl 0.25 --Cl 1.0 --Cq 1.0 --xfilter 0.0 "+\
            "--nPerh 2.01 --graphics False --clearDirectories False --doCompare False "+\
            "--dataDir None --vizTime None --vizCycle None"

    def set_gen_inputs(self):
        steps = 100
        rmin = 0.
        rmax = 1.
        thetaFactor = 0.5
        Ntotal2d = self.ncores*n_per_core_2d
        # Determine nRadial to get Ntotal2d number of SPH nodes for a constantDTheta distribution
        area = np.pi*rmax**2*thetaFactor/2.
        dr = np.sqrt(area/Ntotal2d)
        nRadial = int(rmax/dr)
        self.gen_inp = f"{self.noh_gen_inp} --nRadial {nRadial} --steps {steps}"

#---------------------------------------------------------------------------
# 3D NOH tests
#---------------------------------------------------------------------------
class NOH3D(NOH2D):
    def __init__(self, ncores):
        super().__init__(ncores)
        # Only use half the number of cores
        self.ncores = int(ncores/2)
        self.test_name = "NS3D"
        self.test_file = "functional/Hydro/Noh/Noh-spherical-3d.py"

    def set_gen_inputs(self):
        Ntotal = self.ncores*n_per_core_3d
        npd = int(np.cbrt(Ntotal))
        steps = 10
        self.gen_inp = f"{self.noh_gen_inp} --nx {npd} --ny {npd} --nz {npd} --steps {steps}"

#---------------------------------------------------------------------------
# General functions
#---------------------------------------------------------------------------

# Use arbitrary number of cores if just need test names
def get_all_tests(ncores = 0):
    # Recursively retrieve all subclasses of TestParams
    seen = set()
    work = [TestParams]
    while work:
        parent = work.pop()
        for child in parent.__subclasses__():
            seen.add(child(ncores))
            work.append(child)
    return list(seen)

def get_all_test_names():
    tests = get_all_tests()
    names = []
    for i in tests:
        names.extend(i.test_names())
    return names
