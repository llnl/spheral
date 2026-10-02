"""
Spheral Threading module.  Support for CPU and GPU threading.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Threading/GPUUtils.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#...............................................................................
# init GPUs
# stack_mult is the number of bytes to increase the device stack limit to
@PYB11cppname("GPUUtils::deviceCount")
def deviceCount():
    return "int"

@PYB11cppname("GPUUtils::initGPUs")
def initGPUs(stack_mult = ("int", "8")):
    return "void"

@PYB11cppname("GPUUtils::deviceSync")
def deviceSync():
    return "void"
