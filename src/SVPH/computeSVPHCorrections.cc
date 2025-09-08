//------------------------------------------------------------------------------
// Modified ideas from CRKSPH using the SVPH normalized kernel estimates.
//------------------------------------------------------------------------------
#include "computeSVPHCorrections.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/innerProduct.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Internal utility functions in the unnamed namespace.
//------------------------------------------------------------------------------
namespace {

// 1D
static inline
Dim<1>::ThirdRankTensor
gradxij2(const Dim<1>::Vector& xij) {
  Dim<1>::ThirdRankTensor result;
  result(0,0,0) = 2.0*xij(0);
  return result;
}

// 2D
static inline
Dim<2>::ThirdRankTensor
gradxij2(const Dim<2>::Vector& xij) {
  Dim<2>::ThirdRankTensor result;
  result(0,0,0) = 2.0*xij(0);
  result(0,0,1) =        0.0;
  result(0,1,0) =     xij(1);
  result(0,1,1) =     xij(0);
  result(1,0,0) =     xij(1);
  result(1,0,1) =     xij(0);
  result(1,1,0) =        0.0;
  result(1,1,1) = 2.0*xij(1);
  return result;
}

// 3D
static inline
Dim<3>::ThirdRankTensor
gradxij2(const Dim<3>::Vector& /*xij*/) {
  Dim<3>::ThirdRankTensor result;
  // result(0,0,0) = 2.0*xij(0);
  // result(0,0,1) =     xij(1);
  // result(0,1,0) =     xij(1);
  // result(0,1,1) =        0.0;
  // result(1,0,0) =        0.0;
  // result(1,0,1) =     xij(0);
  // result(1,1,0) =     xij(0);
  // result(1,1,1) = 2.0*xij(1);
  return result;
}

}

//------------------------------------------------------------------------------
// Compute corrections for a single NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeSVPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Scalar>& volume,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       Field<Dimension, typename Dimension::Scalar>& A,
                       Field<Dimension, typename Dimension::Vector>& B,
                       Field<Dimension, typename Dimension::Tensor>& gradB) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Pre-conditions.
  const size_t numNodeLists = volume.size();
  const NodeList<Dimension>& nodeList = A.nodeList();
  const int firstGhostNodei = nodeList.firstGhostNode();
  CONTRACT_VAR(numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(volume.haveNodeList(nodeList));
  REQUIRE(position.haveNodeList(nodeList));
  REQUIRE(H.haveNodeList(nodeList));

  // Zero out the result.
  A = 0.0;
  B = Vector::zero;
  gradB = Tensor::zero;

  // Figure out which NodeList we're working on.
  const size_t nodeListi = std::distance(volume.begin(), volume.fieldForNodeList(nodeList));

  // We can derive everything in terms of the first and second moments 
  // of the local positions.
  Field<Dimension, Vector> m1("first moment", nodeList);
  Field<Dimension, SymTensor> m2("second moment", nodeList);
  Field<Dimension, Tensor> gradm1("gradient of the first moment", nodeList);
  Field<Dimension, ThirdRankTensor> gradm2("gradient of the second moment", nodeList);

  // Iterate over the nodes in this node list.
  for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
       iItr != connectivityMap.end(nodeListi);
       ++iItr) {
    const int i = *iItr;

    // Get the state for node i.
    const Scalar Vi = volume(nodeListi, i);
    const Vector& ri = position(nodeListi, i);
    const SymTensor& Hi = H(nodeListi, i);
    const Scalar Hdeti = Hi.Determinant();

    // Self contribution.
    A(i) += Vi*W(0.0, Hdeti);
    gradm1(i) += Vi*W(0.0, Hdeti) * Tensor::one;

    // Neighbors!
    const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
    CHECK(fullConnectivity.size() == numNodeLists);
    const vector<int>& connectivity = fullConnectivity[nodeListi];

    // Iterate over the neighbors for in this NodeList.
    for (vector<int>::const_iterator jItr = connectivity.begin();
         jItr != connectivity.end();
         ++jItr) {
      const int j = *jItr;

        // // Check if this node pair has already been calculated.
        // if (connectivityMap.calculatePairInteraction(nodeListi, i, 
        //                                              nodeListi, j,
        //                                              firstGhostNodei)) {

      // State of node j.
      const Scalar Vj = volume(nodeListi, j);
      const Vector& rj = position(nodeListi, j);
      const SymTensor& Hj = H(nodeListi, j);
      const Scalar Hdetj = Hj.Determinant();

      // Kernel weighting and gradient.
      const Vector rij = ri - rj;
      const Vector etai = Hi*rij;
      const Vector etaj = Hj*rij;
      Scalar Wi, gWi, Wj, gWj;
      W.kernelAndGradValue(etai.magnitude(), Hdeti, Wi, gWi);
      W.kernelAndGradValue(etaj.magnitude(), Hdetj, Wj, gWj);
      //const Vector gradWi = -(Hi*etai.unitVector())*gWi;
      const Vector gradWj =  (Hj*etaj.unitVector())*gWj;

      // Normalization.
      A(i) += Vj*Wj;

      // First moment. 
      m1(i) += Vj*Wj * rij;
      gradm1(i) += Vj*(Wj*Tensor::one + outerProduct<Dimension>(rij, gradWj));

      // Second moment.
      const SymTensor thpt = rij.selfdyad();
      m2(i) += Vj*Wj * thpt;
      // gradm2(i) += Vj*(outerProduct<Dimension>(gradWj, thpt) +
      //                  Wj*(outerProduct<Dimension>(rij, Tensor::one) + outerProduct<Dimension>(Tensor::one, rij)));
      gradm2(i) += Vj*(Wj*gradxij2(rij) + outerProduct<Dimension>(gradWj, thpt));

          // // First moment. 
          // m1(i) += Vj*Wj * rij;
          // m1(j) -= Vi*Wi * rij;
          // gradm1(i) += Vj*( rij*gradWj + Tensor::one*Wj);
          // gradm1(j) += Vi*(-rij*gradWi + Tensor::one*Wi);

          // // Second moment.
          // const SymTensor thpt = rij.selfdyad();
          // m2(i) += Vj*Wj * thpt;
          // m2(j) += Vi*Wi * thpt;
          // gradm2(i) += Vj*outerProduct<Dimension>(thpt, gradWj);
          // gradm2(j) += Vi*outerProduct<Dimension>(thpt, gradWi);
          // for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
          //   for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
          //     gradm2(i)(ii, jj, jj) += Vj*Wj * rij(ii);
          //     gradm2(j)(ii, jj, jj) -= Vi*Wi * rij(ii);

          //     gradm2(i)(ii, jj, ii) += Vj*Wj * rij(jj);
          //     gradm2(j)(ii, jj, ii) -= Vi*Wi * rij(jj);
          //   }
          // }

        // }
    }

    // Based on the moments we can calculate the SVPH corrections terms and their gradients.
    if (i < firstGhostNodei) {
      CHECK(A(i) > 0.0);
      A(i) = 1.0/A(i);

      CHECK2(abs(m2(i).Determinant()) > 1.0e-30, i << " " << m2(i).Determinant());
      const SymTensor m2inv = m2(i).Inverse();
      B(i) = -(m2inv*m1(i));

      gradB(i) = -innerProduct<Dimension>(m2inv, gradm1(i)) + 
        innerProduct<Dimension>(innerProduct<Dimension>(innerProduct<Dimension>(m2inv, gradm2(i)), m2inv), m1(i));

      // gradB(nodeListi, i) = -innerProduct<Dimension>(gradm1(i), m2inv) + 
      //   innerProduct<Dimension>(innerProduct<Dimension>(innerProduct<Dimension>(m1(i), m2inv), gradm2(i)), m2inv);

      // if (i == 1275) {
      //   cerr << ri << " " << Vi << " " << Hi << endl
      //        << "m1     = " << m1(i) << endl
      //        << "m2     = " << m2(i) << endl
      //        << "gradm1 = " << gradm1(i) << endl
      //        << "gradm2 = " << gradm2(i) << endl
      //        << "B      = " << B(nodeListi, i) << endl
      //        << "gradB  = " << gradB(nodeListi, i) << endl;
      // }

    }
  }
}

//------------------------------------------------------------------------------
// Compute corrections for a all NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeSVPHCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Scalar>& volume,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       FieldList<Dimension, typename Dimension::Scalar>& A,
                       FieldList<Dimension, typename Dimension::Vector>& B,
                       FieldList<Dimension, typename Dimension::Tensor>& gradB) {
  const size_t numNodeLists = A.size();
  for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    computeSVPHCorrections<Dimension>(connectivityMap,
                                      W,
                                      volume,
                                      position,
                                      H,
                                      *A[nodeListi],
                                      *B[nodeListi],
                                      *gradB[nodeListi]);
  }
}

// //------------------------------------------------------------------------------
// // 2D specialization.
// //------------------------------------------------------------------------------
// template<>
// void
// computeSVPHCorrections<Dim<2> >(const ConnectivityMap<Dim<2> >& connectivityMap,
//                                 const TableKernel<Dim<2> >& W,
//                                 const FieldList<Dim<2> , Dim<2>::Scalar>& volume,
//                                 const FieldList<Dim<2> , Dim<2>::Vector>& position,
//                                 const FieldList<Dim<2> , Dim<2>::SymTensor>& H,
//                                 FieldList<Dim<2> , Dim<2>::Scalar>& A,
//                                 FieldList<Dim<2> , Dim<2>::Vector>& B,
//                                 FieldList<Dim<2> , Dim<2>::Tensor>& gradB) {

//   // Pre-conditions.
//   const size_t numNodeLists = A.size();
//   REQUIRE(volume.size() == numNodeLists);
//   REQUIRE(position.size() == numNodeLists);
//   REQUIRE(H.size() == numNodeLists);
//   REQUIRE(B.size() == numNodeLists);
//   REQUIRE(gradB.size() == numNodeLists);

//   typedef Dim<2> Dimension;
//   typedef Dimension::Scalar Scalar;
//   typedef Dimension::Vector Vector;
//   typedef Dimension::Tensor Tensor;
//   typedef Dimension::SymTensor SymTensor;
//   typedef Dimension::ThirdRankTensor ThirdRankTensor;

//   // Zero out the result.
//   A = 0.0;
//   B = Vector::zero;
//   gradB = Tensor::zero;

//   // Walk the FluidNodeLists.
//   for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
//     const NodeList<Dimension>& nodeList = B[nodeListi]->nodeList();
//     const int firstGhostNodei = nodeList.firstGhostNode();

//     // We can derive everything in terms of the first and second moments 
//     // of the local positions.
//     Field<Dimension, Vector> m1("first moment", nodeList);
//     Field<Dimension, SymTensor> m2("second moment", nodeList);
//     Field<Dimension, Tensor> gradm1("gradient of the first moment", nodeList);
//     Field<Dimension, ThirdRankTensor> gradm2("gradient of the second moment", nodeList);
//     Field<Dimension, ThirdRankTensor> phi("phi", nodeList);
//     Field<Dimension, Vector> m2det_1("m2det_1", nodeList);
//     Field<Dimension, Vector> m2det_2("m2det_2", nodeList);
//     Field<Dimension, Vector> m2det_3("m2det_3", nodeList);

//     // Iterate over the nodes in this node list.
//     for (ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
//          iItr != connectivityMap.end(nodeListi);
//          ++iItr) {
//       const int i = *iItr;

//       // Get the state for node i.
//       const Scalar Vi = volume(nodeListi, i);
//       const Vector& ri = position(nodeListi, i);
//       const SymTensor& Hi = H(nodeListi, i);

//       // Self contribution.
//       A(nodeListi, i) += Vi*W(0.0, 1.0);
//       gradm1(i) += Vi*W(0.0, 1.0);

//       // Neighbors!
//       const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
//       CHECK(fullConnectivity.size() == numNodeLists);
//       const vector<int>& connectivity = fullConnectivity[nodeListi];

//       // Iterate over the neighbors for in this NodeList.
//       for (vector<int>::const_iterator jItr = connectivity.begin();
//            jItr != connectivity.end();
//            ++jItr) {
//         const int j = *jItr;

//         // // Check if this node pair has already been calculated.
//         // if (connectivityMap.calculatePairInteraction(nodeListi, i, 
//         //                                              nodeListi, j,
//         //                                              firstGhostNodei)) {

//           // State of node j.
//           const Scalar Vj = volume(nodeListi, j);
//           const Vector& rj = position(nodeListi, j);
//           const SymTensor& Hj = H(nodeListi, j);

//           // Kernel weighting and gradient.
//           const Vector rij = ri - rj;
//           const Vector etai = Hi*rij;
//           const Vector etaj = Hj*rij;
//           const std::pair<double, double> WWi = W.kernelAndGradValue(etai.magnitude(), 1.0);
//           const Scalar& Wi = WWi.first;
//           const Vector gradWi = -(Hi*etai.unitVector())*WWi.second;
//           const std::pair<double, double> WWj = W.kernelAndGradValue(etaj.magnitude(), 1.0);
//           const Scalar& Wj = WWj.first;
//           const Vector gradWj = (Hj*etaj.unitVector())*WWj.second;

//           // Normalization.
//           A(nodeListi, i) += Vj*Wj;

//           // First moment. 
//           m1(i) += Vj*Wj * rij;
//           gradm1(i) += Vj*Wj * Tensor(Wj + rij(0)*gradWj(0),      rij(0)*gradWj(1),
//                                            rij(1)*gradWj(0), Wj + rij(1)*gradWj(1));
//           // gradm1(i) += Vj*(Wj*Tensor::one + outerProduct<Dimension>(rij, gradWj));

//           // Second moment.
//           const SymTensor thpt = rij.selfdyad();
//           m2(i) += Vj*Wj * thpt;
//           gradm2(i)(0,0,0) += Vj*(2.0*rij(0)*Wj + rij(0)*rij(0)*gradWj(0));
//           gradm2(i)(0,0,1) += Vj*(                rij(0)*rij(0)*gradWj(1));
//           gradm2(i)(0,1,0) += Vj*(    rij(1)*Wj + rij(0)*rij(1)*gradWj(0));
//           gradm2(i)(0,1,1) += Vj*(    rij(0)*Wj + rij(0)*rij(1)*gradWj(1));
//           gradm2(i)(1,0,0) += Vj*(    rij(1)*Wj + rij(0)*rij(1)*gradWj(0));
//           gradm2(i)(1,0,1) += Vj*(    rij(0)*Wj + rij(0)*rij(1)*gradWj(1));
//           gradm2(i)(1,1,0) += Vj*(                rij(1)*rij(1)*gradWj(0));
//           gradm2(i)(1,1,1) += Vj*(2.0*rij(1)*Wj + rij(1)*rij(1)*gradWj(1));
//           // gradm2(i) += Vj*(Wj*gradxij2(rij) + outerProduct<Dimension>(gradWj, thpt));

//           phi(i)(0,0,0) += Vj*(                rij(1)*rij(1)*gradWj(0));
//           phi(i)(0,0,1) += Vj*(2.0*rij(1)*Wj + rij(1)*rij(1)*gradWj(1));
//           phi(i)(0,1,0) -= Vj*(    rij(1)*Wj + rij(0)*rij(1)*gradWj(0));
//           phi(i)(0,1,1) -= Vj*(    rij(0)*Wj + rij(0)*rij(1)*gradWj(1));
//           phi(i)(1,0,0) -= Vj*(    rij(1)*Wj + rij(0)*rij(1)*gradWj(0));
//           phi(i)(1,0,1) -= Vj*(    rij(0)*Wj + rij(0)*rij(1)*gradWj(1));
//           phi(i)(1,1,0) += Vj*(2.0*rij(0)*Wj + rij(0)*rij(0)*gradWj(0));
//           phi(i)(1,1,1) += Vj*(                rij(0)*rij(0)*gradWj(1));

//           m2det_1(i) += Vj*Vector(                rij(1)*rij(1)*gradWj(0),
//                                   2.0*Wj*rij(1) + rij(1)*rij(1)*gradWj(1));
//           m2det_2(i) += Vj*Vector(2.0*Wj*rij(0) + rij(0)*rij(0)*gradWj(0),
//                                                   rij(0)*rij(0)*gradWj(1));
//           m2det_3(i) += Vj*Vector(rij(1)*Wj + rij(0)*rij(1)*gradWj(0),
//                                   rij(0)*Wj + rij(0)*rij(1)*gradWj(1));

//         // }
//       }

//       // Based on the moments we can calculate the SVPH corrections terms and their gradients.
//       if (i < firstGhostNodei) {
//         CHECK(A(nodeListi, i) > 0.0);
//         A(nodeListi, i) = 1.0/A(nodeListi, i);

//         CHECK2(abs(m2(i).Determinant()) > 1.0e-30, i << " " << m2(i).Determinant());
//         const SymTensor m2inv = m2(i).Inverse();
//         B(nodeListi, i) = -(m2inv*m1(i));

//         // gradB(nodeListi, i) = -innerProduct<Dimension>(m2inv, gradm1(i)) + 
//         //   innerProduct<Dimension>(innerProduct<Dimension>(innerProduct<Dimension>(m2inv, gradm2(i)), m2inv), m1(i));

//         const Scalar m2det = m2(i).Determinant();
//         const Vector gradm2det = m2(i)(0,0)*m2det_1(i) + m2(i)(1,1)*m2det_2(i) - 2.0*m2(i)(1,0)*m2det_3(i);
//         const SymTensor F( m2(i)(1,1), -m2(i)(0,1),
//                           -m2(i)(1,0),  m2(i)(0,0));
//         const ThirdRankTensor gradm2inv = phi(i)/m2det - outerProduct<Dim<2> >(gradm2det, F)/(m2det*m2det);

//         // gradB(nodeListi, i) = -innerProduct<Dimension>(m2inv, gradm1(i)) + innerProduct<Dimension>(gradm2inv, m1(i));

//         gradB(nodeListi, i)(0,0) = -(F(0,0)*gradm1(i)(0,0) + m1(i)(0)*phi(i)(0,0,0) + F(1,0)*gradm1(i)(1,0) + m1(i)(1)*phi(i)(1,0,0) - (F(0,0)*m1(i)(0) + F(1,0)*m1(i)(1))*gradm2det(0)/m2det)/m2det;
//         gradB(nodeListi, i)(0,1) = -(F(0,0)*gradm1(i)(0,1) + m1(i)(0)*phi(i)(0,0,1) + F(1,0)*gradm1(i)(1,1) + m1(i)(1)*phi(i)(1,0,1) - (F(0,0)*m1(i)(0) + F(1,0)*m1(i)(1))*gradm2det(1)/m2det)/m2det;
//         gradB(nodeListi, i)(1,0) = -(F(1,0)*gradm1(i)(0,0) + m1(i)(0)*phi(i)(1,0,0) + F(1,1)*gradm1(i)(1,0) + m1(i)(1)*phi(i)(1,1,0) - (F(1,0)*m1(i)(0) + F(1,1)*m1(i)(1))*gradm2det(0)/m2det)/m2det;
//         gradB(nodeListi, i)(1,1) = -(F(1,0)*gradm1(i)(0,1) + m1(i)(0)*phi(i)(1,0,1) + F(1,1)*gradm1(i)(1,1) + m1(i)(1)*phi(i)(1,1,1) - (F(1,0)*m1(i)(0) + F(1,1)*m1(i)(1))*gradm2det(1)/m2det)/m2det;

//         if (i == 210) { //1275) {
//           cerr << ri << " " << Vi << " " << Hi << endl
//                << "m1     = " << m1(i) << endl
//                << "m2     = " << m2(i) << endl
//                << "gradm1 = " << gradm1(i) << endl
//                << "gradm2 = " << gradm2(i) << endl
//                << "B      = " << B(nodeListi, i) << endl
//                << "gradB  = " << gradB(nodeListi, i) << endl
//                << "gradB? = " << -innerProduct<Dimension>(m2inv, gradm1(i)) + innerProduct<Dimension>(gradm2inv, m1(i)) << endl;

//           const ThirdRankTensor checkInv = -innerProduct<Dimension>(innerProduct<Dimension>(m2inv, gradm2(i)), m2inv);
//           cerr << "Different ways of measuring gradm2 inverse: " << endl
//                << gradm2inv << endl
//                << checkInv << endl;
//         }

//       }
//     }
//   }
// }

}

