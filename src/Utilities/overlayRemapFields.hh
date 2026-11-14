//------------------------------------------------------------------------------
// Use geometric clipping to remap a set of conserved fields.
// Currently only works single NodeList -> single NodeList, no boundaries.
//------------------------------------------------------------------------------
#ifndef __Spheral_overlayRemapFields__
#define __Spheral_overlayRemapFields__

#include "Boundary/Boundary.hh"
#include "Field/Field.hh"
#include <vector>

namespace Spheral {

template<typename Dimension>
void
overlayRemapFields(const std::vector<Boundary<Dimension>*>& boundaries,
                   const std::vector<Field<Dimension, typename Dimension::Scalar>*>& scalarDonorFields,
                   const std::vector<Field<Dimension, typename Dimension::Vector>*>& vectorDonorFields,
                   const std::vector<Field<Dimension, typename Dimension::Tensor>*>& tensorDonorFields,
                   const std::vector<Field<Dimension, typename Dimension::SymTensor>*>& symTensorDonorFields,
                   std::vector<Field<Dimension, typename Dimension::Scalar>*>& scalarAcceptorFields,
                   std::vector<Field<Dimension, typename Dimension::Vector>*>& vectorAcceptorFields,
                   std::vector<Field<Dimension, typename Dimension::Tensor>*>& tensorAcceptorFields,
                   std::vector<Field<Dimension, typename Dimension::SymTensor>*>& symTensorAcceptorFields,
                   const typename Dimension::Scalar etaMax = 1.0);

}

#endif
