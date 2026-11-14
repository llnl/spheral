from SpheralCompiledPackages import *
import warnings

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SolidNodeList defintion.
#-------------------------------------------------------------------------------
def DEMNodeListFactory(ndim):
    suffix = "{}d".format(ndim)
    DEMNodeList = eval("DEMNodeList" + suffix)
    def factory(name,
                numInternal = 0,
                numGhost = 0,
                hmin = 1.0e-20,
                hmax = 1.0e20,
                hminratio = 0.1,
                nPerh = 1.01,
                maxNumNeighbors = 500,
                neighborSearchBuffer=0.1,
                **unused):
        if unused:
            warnings.warn("makeDEMNodeList{}d ignoring unsupported arguments: {}".format(
                ndim, ", ".join(sorted(unused))), stacklevel=2)
        result = DEMNodeList(name, numInternal, numGhost, 
                             hmin, hmax, hminratio, 
                             nPerh, neighborSearchBuffer, maxNumNeighbors)
        return result
    return factory

#-------------------------------------------------------------------------------
# Create the different SolidNodeLists.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec("makeDEMNodeList{ndim}d = DEMNodeListFactory({ndim})".format(ndim=ndim))
