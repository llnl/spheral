#-------------------------------------------------------------------------------
# ArtificialViscosity pure virtual interface
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11ignore
class ArtificialViscosityAbstractMethods:

    @PYB11const
    def QPiij(self,
              QPiij     = "ReturnType&",
              QPiji     = "ReturnType&",
              Qij       = "Scalar&",
              Qji       = "Scalar&",
              nodeListi = "const unsigned",
              i         = "const unsigned",
              nodeListj = "const unsigned",
              j         = "const unsigned",
              xi        = "const Vector&",
              Hi        = "const SymTensor&",
              etai      = "const Vector&",
              vi        = "const Vector&",
              rhoi      = "const Scalar",
              csi       = "const Scalar",
              xj        = "const Vector&",
              Hj        = "const SymTensor&",
              etaj      = "const Vector&",
              vj        = "const Vector&",
              rhoj      = "const Scalar",
              csj       = "const Scalar",
              fCl       = "const FieldList<%(Dimension)s, Scalar>&",
              fCq       = "const FieldList<%(Dimension)s, Scalar>&",
              DvDx      = "const FieldList<%(Dimension)s, Tensor>&"):
        """All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
Returns the pair values QPiij and QPiji by reference as the first two arguments.
Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
 by the ArtficialViscosity (particularly DvDx)."""
        return "void"
