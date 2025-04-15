from math import *
import mpi

import Spheral
from siloPointmeshDump import *
try:
    from SpheralVoronoiSiloDump import SpheralVoronoiSiloDump
except:
    print("centroidalRelaxNoddes unable to import SpheralVoronoiSiloDump -- no tessellation output supported.")
    SpheralVoronoiSiloDump = None

#-------------------------------------------------------------------------------
# Centroidally (in mass) relax points allowing a linear density gradient.
#-------------------------------------------------------------------------------
def centroidalRelaxNodes(nodeListsAndBounds,
                         W,
                         rho,
                         gradrho = None,
                         boundaries = [],
                         maxIterations = 100,
                         maxFracTol = 1.0e-2,
                         avgFracTol = 1.0e-3,
                         correctionOrder = Spheral.LinearOrder,
                         centroidFrac = 0.25,
                         tessellationBaseDir = ".",
                         tessellationFileName = None):

    # Decide on our dimensionality and import the appropriate aliases.
    assert (isinstance(W, Spheral.TableKernel1d) or
            isinstance(W, Spheral.TableKernel2d) or
            isinstance(W, Spheral.TableKernel3d))
    if isinstance(W, Spheral.TableKernel1d):
        import Spheral1d as sph
        FacetedVolume = sph.Box1d
        ndim = 1
    elif isinstance(W, Spheral.TableKernel2d):
        import Spheral2d as sph
        FacetedVolume = sph.Polygon
        ndim = 2
    else:
        import Spheral3d as sph
        FacetedVolume = sph.Polyhedron
        ndim = 3

    # Did we get passed a function or a constant for the density?
    if type(rho) is float:
        rhoConst = True
        class rhofunctor(sph.VectorScalarFunctor):
            def __init__(self):
                sph.VectorScalarFunctor.__init__(self)
            def __call__(self, posi):
                return rho
    else:
        rhoConst = False
        class rhofunctor(sph.VectorScalarFunctor):
            def __init__(self):
                sph.VectorScalarFunctor.__init__(self)
            def __call__(self, posi):
                return rho(posi)
    rhofunc = rhofunctor()

    # What about the gradrho?  Did we get passed anything?
    if gradrho is None:
        useGradRhoFunc = False
        class gradrhofunctor(sph.VectorVectorFunctor):
            def __init__(self):
                sph.VectorVectorFunctor.__init__(self)
            def __call__(self, posi):
                assert "Hey gradrhofunc unimplemented!"
                return 0.0
    else:
        useGradRhoFunc = True
        if type(gradrho) is float:
            class gradrhofunctor(sph.VectorVectorFunctor):
                def __init__(self):
                    sph.VectorVectorFunctor.__init__(self)
                def __call__(self, posi):
                    return gradrho
        else:
            class gradrhofunctor(sph.VectorVectorFunctor):
                def __init__(self):
                    sph.VectorVectorFunctor.__init__(self)
                def __call__(self, posi):
                    return gradrho(posi)
    gradrhofunc = gradrhofunctor()

    # Split out the NodeLists and bounding volumes (if available), depending on what was passed.
    bounds = sph.vector_of_FacetedVolume()
    holes = sph.vector_of_vector_of_FacetedVolume()
    for x in nodeListsAndBounds:
        if type(nodeListsAndBounds[0]) is tuple:
            bounds = sph.vector_of_FacetedVolume([FacetedVolume()]*len(nodeListsAndBounds))
            holes = sph.vector_of_vector_of_FacetedVolume([sph.vector_of_FacetedVolume()]*len(nodeListsAndBounds))
    if len(bounds) > 0:
        nodeLists = []
        for i, xtup in enumerate(nodeListsAndBounds):
            if type(xtup) is tuple:
                nodeLists.append(xtup[0])
                assert len(xtup) in (2,3)
                bounds[i] = xtup[1]
                if len(xtup) == 3:    # Check for holes
                    assert type(xtup[2]) is list
                    for x in xtup[2]:
                        holes[i].append(x)
            else:
                nodeLists.append(xtup)
    else:
        nodeLists = nodeListsAndBounds

    # Build a local DataBase.
    db = sph.DataBase()
    for nodes in nodeLists:
        db.appendNodeList(nodes)

    # We need the boundaries as a vector
    bound_vec = sph.vector_of_Boundary()
    for bc in boundaries:
        bound_vec.append(bc)

    # Prepare the return FieldLists.
    vol = db.newFluidScalarFieldList(0.0, "volume")
    surfacePoint = sph.IntFieldList()
    cells = sph.FacetedVolumeFieldList()
    cellFaceFlags = db.newGlobalvector_of_CellFaceFlagFieldList(sph.vector_of_CellFaceFlag(), "face flags")
    etaVoidPoints = db.newGlobalvector_of_VectorFieldList(eval("sph.vector_of_Vector%id()" % db.nDim), "eta void points")

    # Initialize volume
    massf = db.fluidMass
    rhof = db.fluidMassDensity
    numNodeLists = db.numFluidNodeLists
    for k in range(numNodeLists):
        n = massf[k].numInternalElements
        for i in range(n):
            assert massf(k,i) > 0.0, "Bad mass (%i,%i), %g" % (k, i, massf(k,i))
            assert rhof(k,i) > 0.0, "Bad density (%i,%i), %g" % (k, i, rhof(k,i))
            vol[k][i] = massf(k,i)/rhof(k,i)

    # We let the C++ method do the heavy lifting.
    iterations = sph.centroidalRelaxNodesImpl(db,
                                              bounds,
                                              holes,
                                              W,
                                              rhofunc,
                                              gradrhofunc,
                                              rhoConst,
                                              useGradRhoFunc,
                                              bound_vec,
                                              maxIterations,
                                              maxFracTol,
                                              avgFracTol,
                                              correctionOrder,
                                              centroidFrac,
                                              vol,
                                              surfacePoint,
                                              cells)

    # Make a final call to computeVoronoiVolume to get the more expensive surfacePoint and cells fields.
    surfacePoint = db.newFluidIntFieldList(0, "surface point")
    deltaMedian = db.newFluidVectorFieldList(sph.Vector.zero, "delta medial position")
    if tessellationFileName:
        cells = db.newFluidFacetedVolumeFieldList(sph.FacetedVolume(), "cells")
    sph.computeVoronoiVolume(db.fluidPosition, 
                             db.fluidHfield, 
                             db.connectivityMap(),
                             sph.SymTensorFieldList(), # no damage
                             bounds,
                             holes,
                             bound_vec,
                             sph.ScalarFieldList(),   # no weights
                             surfacePoint,
                             vol,
                             deltaMedian,
                             etaVoidPoints,
                             cells,
                             cellFaceFlags)

    # Update the masses using rho and volume.
    rho = db.fluidMassDensity
    for k, nodes in enumerate(db.fluidNodeLists):
        n = nodes.numInternalNodes
        mass = nodes.mass()
        for i in range(n):
            assert vol(k,i) > 0.0
            mass[k] = vol(k,i)*rho(k,i)

    # If requested, dump the final info to a diagnostic viz file.
    if tessellationFileName and SpheralVoronoiSiloDump:
        dumper = SpheralVoronoiSiloDump(baseFileName = tessellationFileName,
                                        baseDirectory = tessellationBaseDir,
                                        listOfFieldLists = [vol, surfacePoint, db.fluidMass, db.fluidMassDensity],
                                        boundaries = boundaries,
                                        cells = cells)
        dumper.dump(0.0, iterations)

        siloPointmeshDump(baseName = tessellationFileName + "_points",
                          fieldLists = [vol, surfacePoint, db.fluidMass, db.fluidMassDensity])

    return vol, surfacePoint
