//------------------------------------------------------------------------------
// Compute the CRKSPH corrections.
//------------------------------------------------------------------------------
#include <stdio.h>
#include "computeCRKSPHEvaluation.hh"
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

namespace Spheral {

template<typename Dimension>
void
computeCRKSPHEvaluation(const ConnectivityMap<Dimension>& connectivityMap,
                        const TableKernel<Dimension>& W,
                        const FieldList<Dimension, typename Dimension::Scalar>& weight,
                        const FieldList<Dimension, typename Dimension::Vector>& position,
                        const FieldList<Dimension, typename Dimension::SymTensor>& H,
                        size_t nodeListi, const int i, typename Dimension::Vector reval,
                        const bool coupleNodeLists, typename Dimension::Scalar& WCRKSPH, typename Dimension::Vector& gradWCRKSPH){

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  //YYYY Do I not update ghost nodes? 
  //const int firstGhostNodei = nodeLists[nodeListi]->firstGhostNode();
  const size_t numNodeLists = nodeLists.size();

  // Pre-conditions.
  REQUIRE(weight.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);

  const Scalar wi = weight(nodeListi, i);
  const Vector& ri = position(nodeListi, i);
  const SymTensor& Hi = H(nodeListi, i);
  const Scalar Hdeti = Hi.Determinant();
  const Vector rei = reval - ri;
  const Vector etai = Hi*rei;

  Scalar Wi, gradWSPH;
  W.kernelAndGradValue(etai.magnitude(), Hdeti, Wi, gradWSPH);
  if(abs(WWi.first) < 1.0e-16 and abs(WWi.first) < 1.0e-16){//If SPH kernel and gradient are zero, no point in calculating corrected kernel
    WCRKSPH = 0.0;
    gradWCRKSPH = Vector::zero();
    return;
  }

  // Zero out the result.
  Scalar m0 = 0.0;
  Vector m1 = Vector::zero();
  SymTensor m2 = SymTensor::zero();
  Vector gradm0 = Vector::zero();
  Tensor gradm1 = Tensor::zero();
  ThirdRankTensor gradm2 = ThirdRankTensor::zero();
  Scalar A0 = 0.0;
  Scalar A = 0.0;
  Vector B = Vector::zero();
  Vector C = Vector::zero();
  Tensor D = Tensor::zero();
  Vector gradA0 = Vector::zero();
  Vector gradA = Vector::zero();
  Tensor gradB = Tensor::zero();

  // Neighbors!
  bool first_time=true;//Used as a flag to include self contribution
  const vector<vector<int> > fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
  CHECK(fullConnectivity.size() == numNodeLists);
  for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
    if (coupleNodeLists or nodeListi == nodeListj) {
      const vector<int>& connectivity = fullConnectivity[nodeListj];

      // Iterate over the neighbors for in this NodeList.
      for (vector<int>::const_iterator jItr = connectivity.begin();
           jItr != connectivity.end();
           ++jItr) {
           
          int j = *jItr;
          size_t mynodeList=nodeListj;
          if(first_time){//include self contribution
            j=i;
            --jItr;
            first_time=false;
            mynodeList=nodeListi;
          }
          // State of node j.
          const Scalar wj = weight(mynodeList, j);
          const Vector& rj = position(mynodeList, j);
          const SymTensor& Hj = H(mynodeList, j);
          const Scalar Hdetj = Hj.Determinant();

          // Kernel weighting and gradient.
          const Vector rej = reval - rj;
          const Vector etaj = Hj*rej;
          Scalar Wj, gWj;
          W.kernelAndGradValue(etaj.magnitude(), Hdetj, Wj, gWj);
          const Vector gradWj = (Hj*etaj.unitVector())*gWj;

          // Zeroth moment. 
          const Scalar wwj = wj*Wj;
          m0 += wwj;
          gradm0 += wj*gradWj;

          // First moment. 
          m1 += wwj * rej;
          gradm1 += wj*(outerProduct<Dimension>( rej, gradWj) + Tensor::one()*Wj);

          // Second moment.
          const SymTensor thpt = rej.selfdyad();
          m2 += wwj*thpt;
          gradm2 += wj*outerProduct<Dimension>(thpt, gradWj);

          
          for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
            for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
              gradm2(ii, jj, jj) += wwj*rej(ii);
              gradm2(jj, ii, jj) += wwj*rej(ii);

            }
          }
      }
    }
  }
  // Based on the moments we can calculate the CRKSPH corrections terms and their gradients.
  //if (i < firstGhostNodei) {
  CHECK2(abs(m2.Determinant()) > 1.0e-30, i << " " << m0 << " " << m2 << " " << m2.Determinant());
  const SymTensor m2inv = m2.Inverse();
  const Vector m2invm1 = m2inv*m1;
  const Scalar Ainv = m0 - m2invm1.dot(m1);
  CHECK(Ainv != 0.0);
  A0 = 1.0/m0;
  A = 1.0/Ainv;
  B = -m2invm1;
  gradA0 = -FastMath::square(A0)*gradm0;
  gradA = -A*A*gradm0;
//mjc
/*
  for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
    for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
      for (size_t kk = 0; kk != Dimension::nDim; ++kk) {
         gradA(ii) += A*A*m2inv(jj,kk)*m1(kk)*gradm1(jj,ii);
         gradA(ii) += A*A*m2inv(jj,kk)*gradm1(kk,ii)*m1(jj);
         gradB(ii,jj) -= m2inv(ii,kk)*gradm1(kk,jj);
        for (size_t ll = 0; ll != Dimension::nDim; ++ll) {
          for (size_t mm = 0; mm != Dimension::nDim; ++mm) {
             gradA(ii) -= A*A*m2inv(jj,kk)*gradm2(kk,ll,ii)*m2inv(ll,mm)*m1(mm)*m1(jj);
             gradB(ii,jj) += m2inv(ii,kk)*gradm2(kk,ll,jj)*m2inv(ll,mm)*m1(mm);
          }
        }
      }
    }
  }
*/

  // // BLAGO!
  // // Force only zeroth corrections.
  // A = A0;
  // B = Vector::zero();
  // gradA = -FastMath::square(A)*gradm0;
  // gradB = Tensor::zero();
  // // BLAGO!

  //}

  //Finally Compute Kernel and Gradient 
  WCRKSPH = A*(1.0 + B.dot(rei))*Wi;
  const Vector gradWi = Hi*etai.unitVector() * WWi.second;
  gradWCRKSPH = A*(1.0 + B.dot(rei))*gradWi + A*B*Wi + gradA*(1.0 + B.dot(rei))*Wi;
  for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
    for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
      gradWCRKSPH(ii) += A*Wi*gradB(jj,ii)*rei(jj);
    }
  }


  
}

}
