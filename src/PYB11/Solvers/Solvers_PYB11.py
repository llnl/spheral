"""
Spheral Solvers module.

Provides interfaces to solver libraries wrapped in Spheral.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

# PYB11opaque = ["std::vector<double>"]

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes = ['"Solvers/containsConstantNodes.hh"',
                 '"Solvers/EigenLinearSolver.hh"',
                 '"Solvers/EigenOptions.hh"',
                 '"Solvers/HypreLinearSolver.hh"',
                 '"Solvers/HypreOptions.hh"',
                 '"Solvers/IncrementalStatistic.hh"',
                 '"Solvers/LinearSolver.hh"',
                 '"Solvers/SolverFunction.hh"',
                 '"Solvers/KINSOL.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

from EigenLinearSolver import *
from EigenOptions import *
from HypreLinearSolver import *
from HypreOptions import *
from IncrementalStatistic import *
from KINSOL import *
from LinearSolver import *
from SolverFunction import *

#-------------------------------------------------------------------------------
# Instantiations
#-------------------------------------------------------------------------------
for ndim in dims:
    # Dimension-dependent
    exec('''
@PYB11pycppname("containsConstantNodes")
def containsConstantNodes1%(ndim)id(boundary = "const Boundary<Dim<%(ndim)i>>*"):
    "Does this boundary contain constant boundary nodes?"
    return "bool"
@PYB11pycppname("containsConstantNodes")
def containsConstantNodes2%(ndim)id(boundaries = "const std::vector<Boundary<Dim<%(ndim)i>>*>&"):
    "Do these boundaries contain constant boundary nodes?"
    return "bool"
''' % {"ndim" : ndim})

#-------------------------------------------------------------------------------
# Scalar instantiations
#-------------------------------------------------------------------------------
scalar_types = (("int", "Ordinal"),
                ("double", "Scalar"))

for (value, label) in scalar_types:
    exec('''
%(label)sIncrementalStatistic = PYB11TemplateClass(IncrementalStatistic, template_parameters=("%(value)s"))
''' % {"value" : value,
       "label" : label})
