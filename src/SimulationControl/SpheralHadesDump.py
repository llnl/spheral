#-------------------------------------------------------------------------------
# Dumper class for outputting Spheral++ format data for generating synthetic
# radiographs with Hades.
#-------------------------------------------------------------------------------
import Spheral
from SpheralCompiledPackages import silo
from writeSiloQuadMesh import writeSiloQuadMesh
import mpi
import sys, os, struct, time, bisect
from operator import mul
from functools import reduce

#-------------------------------------------------------------------------------
# Write a silo file resampling to a fixed cartesian mesh for the density.
#-------------------------------------------------------------------------------
def hadesDump(integrator,
              nsample,
              xmin,
              xmax,
              W,
              baseFileName,
              baseDirectory = ".",
              procDirBaseName = "domains",
              mask = None,
              materials = None):

    # Currently suppport 2D and 3D.
    db = integrator.dataBase
    if db.nDim == 2:
        import Spheral2d as sph
    elif db.nDim == 3:
        import Spheral3d as sph
    else:
        raise RuntimeError("hadesDump ERROR: must be 2D or 3D")

    # Prepare to time how long this takes.
    t0 = time.clock()

    # Get the set of material names we're going to write.
    if materials is None:
        materials = list(db.fluidNodeLists)

    # HACK!  We are currently restricting to writing single material output!
    assert len(materials) == 1

    # Make sure the output directory exists.
    if mpi.rank == 0 and not os.path.exists(baseDirectory):
        try:
            os.makedirs(baseDirectory)
        except:
            raise RuntimeError("Cannot create output directory %s" % baseDirectory)
    mpi.barrier()

    # Sample the density.
    ntot = reduce(mul, nsample)
    for nodes in materials:
        print("hadesDump: sampling density for %s..." % nodes.name)
        r = sph.VectorFieldList()
        H = sph.SymTensorFieldList()
        rho = sph.ScalarFieldList()
        r.appendField(nodes.positions())
        H.appendField(nodes.Hfield())
        rho.appendField(nodes.massDensity())

        mf = nodes.mass()
        rhof = nodes.massDensity()
        wf = sph.ScalarField("volume", nodes)
        for i in range(nodes.numNodes):
            wf[i] = mf[i]/max(1e-100, rhof[i])
        w = sph.ScalarFieldList()
        w.copyFields()
        w.appendField(wf)
        #w.appendField(sph.ScalarField("weight", nodes, 1.0))

        fieldListSet = sph.FieldListSet()
        fieldListSet.ScalarFieldLists.append(rho)
        localMask = sph.IntFieldList()
        if mask is None:
            localMask.copyFields()
            localMask.appendField(sph.IntField("mask", nodes, 1))
        else:
            localMask.appendField(mask.fieldForNodeList(nodes))

        scalar_samples = sph.vector_of_vector_of_double()
        vector_samples = sph.vector_of_vector_of_Vector()
        tensor_samples = sph.vector_of_vector_of_Tensor()
        symTensor_samples = sph.vector_of_vector_of_SymTensor()

        (scalar_samples,
         vector_samples,
         tensor_samples,
         symTensor_samples) = sph.sampleMultipleFields2Lattice(fieldListSet,
                                                               r, w, H, localMask,
                                                               W,
                                                               xmin, xmax,
                                                               sph.vector_of_int(nsample))
        print("Generated %i scalar fields" % len(scalar_samples))

        # Write out the silo info
        writeSiloQuadMesh(scalar_data = scalar_samples,
                          ndim = db.nDim,
                          xmin = xmin,
                          xmax = xmax,
                          nglobal = nsample,
                          filename = baseFileName,
                          dirname = baseDirectory,
                          scalar_names = ("den",),
                          materials = materials,
                          time = integrator.currentTime,
                          cycle = integrator.currentCycle,
                          RZ = (GeometryRegistrar.coords() == CoordinateType.RZ))

    mpi.barrier()
    print("hadesDump finished: required %0.2f seconds" % (time.clock() - t0))
    return
