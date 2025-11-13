from Spheral import *
from SpheralTestUtilities import fuzzyEqual
import mpi

#===============================================================================
# Base test class -- generic for all dimensions.
#===============================================================================
class TestSampleMultipleFields2Lattice:

    #---------------------------------------------------------------------------
    # Run the test!
    #---------------------------------------------------------------------------
    def testSample(self):

        if self.ndim == 1:
            from Spheral1d import (ScalarFieldList, VectorFieldList, TensorFieldList, SymTensorFieldList, 
                                   FieldListSet, sampleMultipleFields2LatticeMash,
                                   vector_of_vector_of_Vector, vector_of_vector_of_Tensor, vector_of_vector_of_SymTensor)
        elif self.ndim == 2:
            from Spheral2d import (ScalarFieldList, VectorFieldList, TensorFieldList, SymTensorFieldList, 
                                   FieldListSet, sampleMultipleFields2LatticeMash,
                                   vector_of_vector_of_Vector, vector_of_vector_of_Tensor, vector_of_vector_of_SymTensor)
        else:
            assert self.ndim == 3
            from Spheral3d import (ScalarFieldList, VectorFieldList, TensorFieldList, SymTensorFieldList, 
                                   FieldListSet, sampleMultipleFields2LatticeMash,
                                   vector_of_vector_of_Vector, vector_of_vector_of_Tensor, vector_of_vector_of_SymTensor)

        # Build the FieldLists we're going to sample.
        rho = ScalarFieldList()
        eps = ScalarFieldList()
        vel = VectorFieldList()
        Hfl = SymTensorFieldList()
        rho.appendField(self.nodes.massDensity())
        eps.appendField(self.nodes.specificThermalEnergy())
        vel.appendField(self.nodes.velocity())
        Hfl.appendField(self.nodes.Hfield())
        gradv = self.db.newFluidTensorFieldList(self.gradv0, "fake vel gradient")

        # Put them into a FieldListSet.
        fieldListSet = FieldListSet()
        fieldListSet.ScalarFieldLists.append(rho)
        fieldListSet.ScalarFieldLists.append(eps)
        fieldListSet.VectorFieldLists.append(vel)
        fieldListSet.SymTensorFieldLists.append(Hfl)
        fieldListSet.TensorFieldLists.append(gradv)
        assert len(fieldListSet.ScalarFieldLists) == 2
        assert len(fieldListSet.VectorFieldLists) == 1
        assert len(fieldListSet.SymTensorFieldLists) == 1
        assert len(fieldListSet.TensorFieldLists) == 1

        # Build the mask.
        mask = self.db.newGlobalIntFieldList(1)

        # Now sample those suckers.
        r = VectorFieldList()
        w = ScalarFieldList()
        H = SymTensorFieldList()
        r.appendField(self.nodes.positions())
        w.appendField(self.nodes.mass())
        H.appendField(self.nodes.Hfield())
        scalar_samples, vector_samples, tensor_samples, symtensor_samples = sampleMultipleFields2LatticeMash(fieldListSet,
                                                                                                             r,
                                                                                                             w,
                                                                                                             H,
                                                                                                             mask,
                                                                                                             self.WT,
                                                                                                             self.xmin,
                                                                                                             self.xmax,
                                                                                                             self.nsample)

        # Did we get back the correct numbers of sampled values?
        assert len(scalar_samples) == 2
        assert len(vector_samples) == 1
        assert len(tensor_samples) == 1
        assert len(symtensor_samples) == 1
        ntot = 1
        for x in self.nsample:
            ntot *= x
        for array in scalar_samples:
            assert mpi.allreduce(len(array), mpi.SUM) == ntot
        for array in vector_samples:
            assert mpi.allreduce(len(array), mpi.SUM) == ntot
        for array in tensor_samples:
            assert mpi.allreduce(len(array), mpi.SUM) == ntot
        for array in symtensor_samples:
            assert mpi.allreduce(len(array), mpi.SUM) == ntot
        
        # Extract the arrays.
        rhoarray = scalar_samples[0]
        epsarray = scalar_samples[1]
        velarray = vector_samples[0]
        Harray = symtensor_samples[0]
        gradvelarray = tensor_samples[0]

        # See if the resulting fields are constant enough.
        self.assertTrue(min([fuzzyEqual(x, self.rho0, self.scalarTol) for x in rhoarray]),
                        "Failing rho comparison: expect=%g min=%g, max=%g, tol=%g" % (self.rho0,
                                                                                      min(rhoarray),
                                                                                      max(rhoarray),
                                                                                      self.scalarTol))
        self.assertTrue(min([fuzzyEqual(x, self.eps0, self.scalarTol) for x in epsarray]),
                        "Failing eps comparison: expect=%g min=%g, max=%g, tol=%g" % (self.eps0,
                                                                                      min(epsarray),
                                                                                      max(epsarray),
                                                                                      self.scalarTol))
        self.assertTrue(min([fuzzyEqual((x - self.v0).magnitude(), 0.0, self.vectorTol) for x in velarray]),
                        "Failing vel comparison: min=%g, max=%g, tol=%g" % (min((x - self.v0).magnitude() for x in velarray),
                                                                            max((x - self.v0).magnitude() for x in velarray),
                                                                            self.vectorTol))
        self.assertTrue(min([fuzzyEqual((x - self.gradv0).doubledot(x - self.gradv0), 0.0, self.tensorTol) for x in gradvelarray]),
                        "Failing gradv comparison: min=%g, max=%g, tol=%g" % (min([(x - self.gradv0).doubledot(x - self.gradv0) for x in gradvelarray]),
                                                                              max([(x - self.gradv0).doubledot(x - self.gradv0) for x in gradvelarray]),
                                                                              self.tensorTol))
        self.assertTrue(min([fuzzyEqual((x - self.H0).doubledot(x - self.H0), 0.0, self.tensorTol) for x in Harray]),
                        "Failing H comparison: min=%g, max=%g, tol=%g" % (min([(x - self.H0).doubledot(x - self.H0) for x in Harray]),
                                                                          max([(x - self.H0).doubledot(x - self.H0) for x in Harray]),
                                                                          self.tensorTol))

        return

