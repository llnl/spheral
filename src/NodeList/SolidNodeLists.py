from SpheralCompiledPackages import *
import warnings

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# The generic SolidNodeList defintion.
#-------------------------------------------------------------------------------
def SolidNodeListFactory(ndim):
    suffix = "{}d".format(ndim)
    SolidNodeList = eval("SolidNodeList" + suffix)
    NullStrength = eval("NullStrength" + suffix)
    def factory(name,
                eos,
                strength = NullStrength(),
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
            warnings.warn("makeSolidNodeList{}d ignoring unsupported arguments: {}".format(
                ndim, ", ".join(sorted(unused))), stacklevel=2)
        result = SolidNodeList(name = name,
                               eos = eos,
                               strength = strength,
                               numInternal = numInternal,
                               numGhost = numGhost,
                               hmin = hmin,
                               hmax = hmax,
                               hminratio = hminratio, 
                               nPerh = nPerh,
                               maxNumNeighbors = maxNumNeighbors,
                               rhoMin = rhoMin,
                               rhoMax = rhoMax)
        result.eos = eos
        result.strength = strength
        return result
    return factory

#-------------------------------------------------------------------------------
# Create the different SolidNodeLists.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec("makeSolidNodeList{ndim}d = SolidNodeListFactory({ndim})".format(ndim=ndim))
