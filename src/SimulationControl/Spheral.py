# Master import file to import the Spheral packages and a standard set of
# helper extensions.

# Modified version to be compatible with the pybindgen version of Spheral++.

# ------------------------------------------------------------------------------
# Load up MPI.
# ------------------------------------------------------------------------------
import mpi

import SpheralConfigs
from SpheralUtilities import BuildData
if SpheralConfigs.enable_shared() and \
   not BuildData.cxx_compiler_id == "GNU":
    try:
        import sys, ctypes
        sys.setdlopenflags(sys.getdlopenflags() | ctypes.RTLD_GLOBAL)
    except:
        print("WARNING: unable to set python dl flags on Spheral import.")
        pass

# ------------------------------------------------------------------------------
# Import a scipy module to initialize scipy's shared qhull library before
# spheral's static qhull library.
# ------------------------------------------------------------------------------
import scipy.spatial

# ------------------------------------------------------------------------------
# Import the compiled packages.
# ------------------------------------------------------------------------------
from SpheralCompiledPackages import *

# Aliases for a few objects
FacetedVolume1d = Box1d
FacetedVolume2d = Polygon
FacetedVolume3d = Polyhedron

# ------------------------------------------------------------------------------
# Import the Material python extensions.
# ------------------------------------------------------------------------------
from MaterialUnits import *
from MaterialEquationsOfState import *

# ------------------------------------------------------------------------------
# Import the various NodeLists.
# ------------------------------------------------------------------------------
from FluidNodeLists import *
from SolidNodeLists import *
from VoidNodeLists import *
from DEMNodeLists import *

# ------------------------------------------------------------------------------
# Import hydro configurations
# ------------------------------------------------------------------------------
hydroImports = SpheralConfigs.hydro_imports()
for x in hydroImports:
    exec(f"from {x} import *")

# ------------------------------------------------------------------------------
# Import the SolidMaterial python extensions.
# ------------------------------------------------------------------------------
from SolidMaterialUnits import *
from SolidMaterialEquationsOfState import *

from GradyKippTensorDamage import *
from JohnsonCookDamageFactories import (JohnsonCookDamageConstant,
                                        JohnsonCookDamageGaussian,
                                        JohnsonCookDamageWeibull)

# These are for backwards compatibility and should be deprecated
from PorousEquationOfState import *
from PorousStrengthModel import *

# ------------------------------------------------------------------------------
# Helpful things with strings.
# ------------------------------------------------------------------------------
## import CXXTypes
## def add_string_str(self, pystr):
##     return str(self).join(pystr)

## def add_str_string(self, cppstr):
##     return self.join(str(cppstr))

## CXXTypes.string.__add__ = add_string_str

# ------------------------------------------------------------------------------
# Import the picking helpers that register how to pickle Spheral++ objects.
# ------------------------------------------------------------------------------
from SpheralPickle import *

# ------------------------------------------------------------------------------
# Import the python enhanced FlatFileIO object.
# ------------------------------------------------------------------------------
#from ExtendFlatFileIO import *

# ------------------------------------------------------------------------------
# Import the controller and a standard timer class.
# ------------------------------------------------------------------------------
from SpheralTimer import *
from SpheralController import *

# ------------------------------------------------------------------------------
# Import the command line interpreter.
# ------------------------------------------------------------------------------
from SpheralOptionParser import commandLine

# ------------------------------------------------------------------------------
# See if we can import the polytope bindings.
# ------------------------------------------------------------------------------
try:
    import polytope
except:
    print("WARNING: unable to import polytope python bindings.")

# ------------------------------------------------------------------------------
# Import our shadow layers for augmenting C++ types.
# ------------------------------------------------------------------------------
for shadowedthing in ("TillotsonEquationOfState",
                      "GruneisenEquationOfState",
                      "ConstantStrength",
                      "ProbabilisticDamageModel",
                      "IvanoviSALEDamageModel",
                      "StrainPorosity",
                      "PalphaPorosity",
                      "ANEOS"):
    for dim in dims:
        exec(f"from Shadow{shadowedthing} import {shadowedthing}{dim}d")

# ------------------------------------------------------------------------------
# Output some useful Spheral configuration info to stdout
# ------------------------------------------------------------------------------
print("/------------------------------------------------------------------------------\\")
import Spheral_banner
print("|  %-76s|" % ("  number of MPI tasks       : " + str(mpi.procs)))
print("|  %-76s|" % ("  number of threads per rank: " + str(omp_get_num_threads())))
print("\\------------------------------------------------------------------------------/")

# ------------------------------------------------------------------------------
# Set the prompt just to clear to folks they now have Spheral.
# To maintain sanity by default only have one process print the prompt...
# ------------------------------------------------------------------------------
if mpi.rank == 0:
    sys.ps1 = "Spheral> "
else:
    sys.ps1 = ""

    
