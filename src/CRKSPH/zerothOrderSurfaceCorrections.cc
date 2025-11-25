//------------------------------------------------------------------------------
// Look for any points that touch a surface (multi-material or void).
// For such points:
//   - enforce only zeroth order corrections.
//------------------------------------------------------------------------------
#include "zerothOrderSurfaceCorrections.hh"

namespace Spheral {

template<typename Dimension>
void
zerothOrderSurfaceCorrections(FieldList<Dimension, typename Dimension::Scalar>& A,
                              FieldList<Dimension, typename Dimension::Vector>& B,
                              FieldList<Dimension, typename Dimension::Tensor>& C,
                              FieldList<Dimension, typename Dimension::Vector>& gradA,
                              FieldList<Dimension, typename Dimension::Tensor>& gradB,
                              FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC,
                              const FieldList<Dimension, typename Dimension::Scalar>& m0,
                              const FieldList<Dimension, typename Dimension::Vector>& gradm0,
                              const FieldList<Dimension, int>& surfacePoint) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  const auto numNodeLists = surfacePoint.numFields();
  const auto fixC = C.numFields() == numNodeLists;
  REQUIRE(A.numFields() == numNodeLists);
  REQUIRE(B.numFields() == numNodeLists);
  REQUIRE(gradA.numFields() == numNodeLists);
  REQUIRE(gradB.numFields() == numNodeLists);
  REQUIRE(m0.numFields() == numNodeLists);
  REQUIRE(gradm0.numFields() == numNodeLists);

  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = surfacePoint[nodeListi]->numInternalElements();

#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      if (surfacePoint(nodeListi, i) != 0) {
        CHECK(m0(nodeListi, i) != 0.0);
        A(nodeListi, i) = 1.0/m0(nodeListi, i);
        B(nodeListi, i) = Vector::zero();
        if (fixC) C(nodeListi, i) = Tensor::zero();
        gradA(nodeListi, i) = -FastMath::square(A(nodeListi, i))*gradm0(nodeListi, i);
        gradB(nodeListi, i) = Tensor::zero();
        if (fixC) gradC(nodeListi, i) = ThirdRankTensor::zero();
      }
    }
  }
}

}
