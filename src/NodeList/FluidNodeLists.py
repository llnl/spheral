from SpheralCompiledPackages import *
import warnings

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic FluidNodeList pattern.
#-------------------------------------------------------------------------------
def FluidNodeListFactory(ndim):
    suffix = "{}d".format(ndim)
    FluidNodeList = eval("FluidNodeList" + suffix)
    def factory(name,
                eos,
                numInternal = 0,
                numGhost = 0,
                hmin = 1.0e-20,
                hmax = 1.0e20,
                hminratio = 0.1,
                nPerh = 2.01,
                maxNumNeighbors = 500,
                rhoMin = 1.0e-10,
                rhoMax = 1e10,
                **unused):
        if unused:
            warnings.warn("makeFluidNodeList{}d ignoring unsupported arguments: {}".format(
                ndim, ", ".join(sorted(unused))), stacklevel=2)
        result = FluidNodeList(name, eos, numInternal, numGhost, 
                               hmin, hmax, hminratio, 
                               nPerh, maxNumNeighbors,
                               rhoMin, rhoMax)
        result.eos = eos
        return result
    return factory

#-------------------------------------------------------------------------------
# Create the different dimension implementations.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec("makeFluidNodeList{ndim}d = FluidNodeListFactory({ndim})".format(ndim=ndim))
