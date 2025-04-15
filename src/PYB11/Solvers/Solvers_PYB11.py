"""
Spheral Solvers module.

Provides interfaces to solver libraries wrapped in Spheral.
"""

from PYB11Generator import *

PYB11opaque = ["std::vector<double>"]

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes = ['"Solvers/SolverFunction.hh"',
                 '"Solvers/KINSOL.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

from SolverFunction import *
from KINSOL import *
