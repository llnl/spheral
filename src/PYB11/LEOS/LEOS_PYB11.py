"""
Spheral LEOS module.

Provides access to the Livermore EOS package.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#from LEOS_bundle import *
from LEOS import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/GeomPlane.hh"',
                  '"Field/Field.hh"',
                  '"Material/PhysicalConstants.hh"',
                  '"Material/EquationOfState.hh"',
                  #'"LEOS/LEOS_bundle.hh"',
                  '"LEOS/LEOS.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
@PYB11implementation('[](const std::vector<int> eosNumbers, const std::vector<std::string> libNames) { std::cerr << "DEPRECATION WARNING: LEOS_initialize no longer required (no-op)" << std::endl; }')
def LEOS_initialize(eosNumbers = "const std::vector<int>",
                    libNames = ("const std::vector<std::string>", "std::vector<std::string>()")):
    "DEPRECATED -- no longer required"
    return "void"

@PYB11implementation('[](py::list eosNumbers, py::list libNames) { std::cerr << "DEPRECATION WARNING: LEOS_initialize no longer required (no-op)" << std::endl; }')
def LEOS_initialize(eosNumbers = "py::list",
                    libNames = ("py::list", "py::list()")):
    "DEPRECATED -- no longer required"
    return "void"

#-------------------------------------------------------------------------------
# Dimensional instantiations
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
LEOS%(ndim)id = PYB11TemplateClass(LEOS, template_parameters="%(Dimension)s")
''' % {"ndim" : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
