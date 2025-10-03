text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Utilities/overlayRemapFields.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
template void 
overlayRemapFields(const std::vector<Boundary<Dim<%(ndim)s>>*>& boundaries,
                   const std::vector<Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>*>& scalarDonorFields,
                   const std::vector<Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>*>& vectorDonorFields,
                   const std::vector<Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>*>& tensorDonorFields,
                   const std::vector<Field<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>*>& symTensorDonorFields,
                   std::vector<Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Scalar>*>& scalarAcceptorFields,
                   std::vector<Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Vector>*>& vectorAcceptorFields,
                   std::vector<Field<Dim<%(ndim)s>, Dim<%(ndim)s>::Tensor>*>& tensorAcceptorFields,
                   std::vector<Field<Dim<%(ndim)s>, Dim<%(ndim)s>::SymTensor>*>& symTensorAcceptorFields);
}

"""

# Only supports 2D and 3D
dimensions = [2, 3]
