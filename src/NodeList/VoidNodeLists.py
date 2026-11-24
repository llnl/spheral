from SpheralCompiledPackages import *
import warnings

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic VoidNodeList pattern.
#-------------------------------------------------------------------------------
def VoidNodeListFactory(ndim):
    suffix = "{}d".format(ndim)
    NodeList = eval("NodeList" + suffix)
    def factory(name,
                numInternal = 0,
                numGhost = 0,
                hmin = 1.0e-20,
                hmax = 1.0e20,
                hminratio = 0.1,
                nPerh = 2.01,
                maxNumNeighbors = 500,
                **unused):
        if unused:
            warnings.warn("makeVoidNodeList{}d ignoring unsupported arguments: {}".format(
                ndim, ", ".join(sorted(unused))), stacklevel=2)
        result = NodeList(name, numInternal, numGhost, 
                          hmin, hmax, hminratio, 
                          nPerh, maxNumNeighbors)
        return result
    return factory

#-------------------------------------------------------------------------------
# Create the different dimension implementations.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec("makeVoidNodeList{ndim}d = VoidNodeListFactory({ndim})".format(ndim=ndim))
