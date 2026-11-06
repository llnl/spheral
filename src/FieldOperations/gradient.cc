//---------------------------------Spheral++----------------------------------//
// FieldListFunctions -- A set of global functions which can be applied to
// FieldLists.
//
// Created by JMO, Wed Dec  6 21:09:29 PST 2000
//----------------------------------------------------------------------------//
#include "FieldListFunctions.hh"
#include "FieldOperations/FieldListConnectivity.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "Kernel/TableKernel.hh"
#include "Geometry/MathTraits.hh"

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
// Calculate the gradient of a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradient(const FieldList<Dimension, DataType>& fieldList,
         const FieldList<Dimension, typename Dimension::Vector>& position,
         const FieldList<Dimension, typename Dimension::Scalar>& /*weight*/,
         const FieldList<Dimension, typename Dimension::Scalar>& mass,
         const FieldList<Dimension, typename Dimension::Scalar>& rho,
         const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
         const TableKernel<Dimension>& kernel) {

  // Typedef's to ease typing/understandability.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Return FieldList.
  FieldList<Dimension, GradientType> result;
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, GradientType>("grad", (*fieldItr)->nodeList()));
  }
  FieldListConnectivity<Dimension, DataType> connectivity(fieldList, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const Vector& ri = position(fieldi, i);
      const SymTensor& Hi = Hfield(fieldi, i);
      const Scalar& rhoi = rho(fieldi, i);
      const DataType& fi = fieldList(fieldi, i);
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
          const Vector& rj = position(fieldj, j);
          const SymTensor& Hj = Hfield(fieldj, j);
          const Scalar& mj = mass(fieldj, j);
          const DataType& fj = fieldList(fieldj, j);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Vector etaiNorm = etai.unitVector();
          const Vector etajNorm = etaj.unitVector();

          // Get the symmetrized kernel gradient for this node pair.
          Vector gradWij;
          switch(searchType) {
          case NeighborSearchType::GatherScatter:
            gradWij = 0.5*(Hi*etaiNorm*kernel.grad(etai.magnitude(), Hi.Determinant()) + 
                           Hj*etajNorm*kernel.grad(etaj.magnitude(), Hj.Determinant()));
            break;

          case NeighborSearchType::Gather:
            gradWij = Hi*etaiNorm*kernel.grad(etai.magnitude(), Hi.Determinant());
            break;

          case NeighborSearchType::Scatter:
            gradWij = Hj*etajNorm*kernel.grad(etaj.magnitude(), Hj.Determinant());
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Add this nodes contribution to the master value.
          result(fieldi, i) += mj*(fj - fi)*gradWij;
      });

        // Normalize by density.
        CHECK(rhoi > 0.0);
      result(fieldi, i) /= rhoi;
    }
  }

  return result;
}

template<typename Dimension, typename DataType>
FieldList<Dimension, std::vector<typename MathTraits<Dimension, DataType>::GradientType>>
gradient(const FieldList<Dimension, std::vector<DataType>>& fieldList,
         const FieldList<Dimension, typename Dimension::Vector>& position,
         const FieldList<Dimension, typename Dimension::Scalar>& /*weight*/,
         const FieldList<Dimension, typename Dimension::Scalar>& mass,
         const FieldList<Dimension, typename Dimension::Scalar>& rho,
         const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
         const TableKernel<Dimension>& kernel) {

  // Typedef's to ease typing/understandability.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Get size of vector
  size_t vectorSize = 0u;
  bool foundInternalNode = false;
  for (auto fieldi = 0u; fieldi < fieldList.numFields() and not foundInternalNode; ++fieldi) {
    if (fieldList[fieldi]->nodeListPtr()->numInternalNodes() > 0) {
      vectorSize = fieldList(fieldi, 0).size();
      foundInternalNode = true;
    }
  }
  
  // Return FieldList.
  FieldList<Dimension, std::vector<GradientType>> result;
  result.copyFields();
  for (typename FieldList<Dimension, std::vector<DataType>>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, std::vector<GradientType>>("grad",
                                                                   (*fieldItr)->nodeList(),
                                                                   std::vector<GradientType>(vectorSize, GradientType::zero())));
  }
  FieldListConnectivity<Dimension, std::vector<DataType>> connectivity(fieldList, kernel.kernelExtent());
  const auto searchType = connectivity.searchType();
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
      const Vector& ri = position(fieldi, i);
      const SymTensor& Hi = Hfield(fieldi, i);
      const Scalar& rhoi = rho(fieldi, i);
      const auto& fi = fieldList(fieldi, i);
      connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
          const Vector& rj = position(fieldj, j);
          const SymTensor& Hj = Hfield(fieldj, j);
          const Scalar& mj = mass(fieldj, j);
          const auto& fj = fieldList(fieldj, j);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Vector etaiNorm = etai.unitVector();
          const Vector etajNorm = etaj.unitVector();

          // Get the symmetrized kernel gradient for this node pair.
          Vector gradWij;
          switch(searchType) {
          case NeighborSearchType::GatherScatter:
            gradWij = 0.5*(Hi*etaiNorm*kernel.grad(etai.magnitude(), Hi.Determinant()) + 
                           Hj*etajNorm*kernel.grad(etaj.magnitude(), Hj.Determinant()));
            break;

          case NeighborSearchType::Gather:
            gradWij = Hi*etaiNorm*kernel.grad(etai.magnitude(), Hi.Determinant());
            break;

          case NeighborSearchType::Scatter:
            gradWij = Hj*etajNorm*kernel.grad(etaj.magnitude(), Hj.Determinant());
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Add this nodes contribution to the master value.
          for (auto m = 0u; m < vectorSize; ++m) {
            CHECK(result(fieldi, i).size() == vectorSize &&
                  fj.size() == vectorSize &&
                  fi.size() == vectorSize);
            result(fieldi, i)[m] += mj*(fj[m] - fi[m])*gradWij;
          }
      });

        // Normalize by density.
        CHECK(rhoi > 0.0);
        for (auto m = 0u; m < vectorSize; ++m) {
          CHECK(result(fieldi, i).size() == vectorSize);
          result(fieldi, i)[m] /= rhoi;
        }
    }
  }

  return result;
}

//------------------------------------------------------------------------------
// Calculate the gradient of a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
void
gradientPairs(FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>& result,
              const FieldList<Dimension, DataType>& field,
              const FieldList<Dimension, typename Dimension::Vector>& position,
              const FieldList<Dimension, typename Dimension::Scalar>& weight,
              const FieldList<Dimension, typename Dimension::SymTensor>& H,
              const ConnectivityMap<Dimension>& conn,
              const TableKernel<Dimension>& kernel) {
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;
  result = GradientType::zero();
  
  const auto& pairs = conn.nodePairList();
  const auto  npairs = pairs.size();

  for (auto k = 0u; k < npairs; ++k) {
      const auto ni = pairs[k].i_list;
      const auto nj = pairs[k].j_list;
      const auto i = pairs[k].i_node;
      const auto j = pairs[k].j_node;
      
      const auto& ri = position(ni, i);
      const auto& rj = position(nj, j);
      const auto& vi = weight(ni, i);
      const auto& vj = weight(nj, j);
      const auto& hi = H(ni, i);
      const auto& hj = H(nj, j);
      const auto& fi = field(ni, i);
      const auto& fj = field(nj, j);

      const auto rij = ri - rj;
      const auto hdeti = hi.Determinant();
      const auto hdetj = hj.Determinant();
      const auto etai = hi * rij;
      const auto etaj = hj * rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      const auto etaUniti = etai.unitVector();
      const auto etaUnitj = etaj.unitVector();
      const auto hetaUniti = hi * etaUniti;
      const auto hetaUnitj = hj * etaUnitj;
      
      const auto dwi = hetaUniti * kernel.gradValue(etaMagi, hdeti);
      const auto dwj = hetaUnitj * kernel.gradValue(etaMagj, hdetj);
      const auto dwij = 0.5 * (dwi + dwj);

      result(ni, i) += vj * (fj - fi) * dwij;
      result(nj, j) -= vi * (fi - fj) * dwij;
  }
}

//------------------------------------------------------------------------------
// The limiter method.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
monotonicLimiter(const typename Dimension::Scalar& dF, 
                 const typename Dimension::Scalar& dFproj,
                 const typename Dimension::Vector& nhat,
                 const double fuzz = 1.0e-15) {
  CONTRACT_VAR(nhat);
  REQUIRE(fuzz > 0.0);
  REQUIRE(fuzzyEqual(nhat.magnitude2(), 1.0));
  const double dFproj2 = dFproj*dFproj;
  return max(0.0, min(1.0, max(dF*dFproj/(dFproj2 + fuzz), fuzz/(dFproj2 + fuzz))));
}

template<typename Dimension>
inline
double
monotonicLimiter(const typename Dimension::Vector& dF,
                 const typename Dimension::Vector& dFproj,
                 const typename Dimension::Vector& nhat,
                 const double fuzz = 1.0e-15) {
  REQUIRE(fuzz > 0.0);
  REQUIRE(fuzzyEqual(nhat.magnitude2(), 1.0));
  return monotonicLimiter<Dimension>(dF.dot(nhat), dFproj.dot(nhat), nhat, fuzz);
}

//------------------------------------------------------------------------------
// Apply a monotonic limiter to a pre-computed gradient.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
FieldList<Dimension, typename Dimension::SymTensor>
limiter(const FieldList<Dimension, DataType>& fieldList,
        const FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>& gradient,
        const FieldList<Dimension, typename Dimension::Vector>& position,
        const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
        const TableKernel<Dimension>& kernel) {

  const int maxIterations = 10;
  const double W0 = kernel(0.0, 1.0);
  const double tiny = 1.0e-15;

  // Typedef's to ease typing/understandability.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Return FieldList.
  FieldList<Dimension, SymTensor> result;
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, SymTensor>("limiter", (*fieldItr)->nodeList()));
  }
  FieldListConnectivity<Dimension, DataType> connectivity(fieldList, kernel.kernelExtent());
  for (auto fieldi = 0u; fieldi < fieldList.numFields(); ++fieldi) {
    const auto numInternal = fieldList[fieldi]->nodeListPtr()->numInternalNodes();
    for (auto i = 0u; i < numInternal; ++i) {
        const Vector& ri = position(fieldi, i);
        const SymTensor& Hi = Hfield(fieldi, i);
        const DataType& fi = fieldList(fieldi, i);
        const GradientType& gradi = gradient(fieldi, i);

        // Prepare the result for this node.
        SymTensor phi = SymTensor::one();

        // We iterate on the limiter.
        int iter = 0;
        Scalar phimin = 0.0;
        while (not fuzzyEqual(phimin, 1.0, 1.0e-5) and iter < maxIterations) {

          // Loop over the refined neighbors.
          Scalar weightSum = 0.0;
          SymTensor phii;
          phimin = 1.0;
          connectivity.forEachNeighbor(fieldi, i, [&](const size_t fieldj, const size_t j) {
            const Vector& rj = position(fieldj, j);
            const DataType& fj = fieldList(fieldj, j);

            // Compute the pair-wise limiting needed.
            const Vector rji = rj - ri;
            const DataType dF = fj - fi;
            const DataType dFproj = (phi*gradi).dot(rji);
            const Vector rhat = rji.unitVector();
            const Scalar phiij = monotonicLimiter<Dimension>(dF, dFproj, rhat);

            // Increment this iterations estimate of the limiter.
            const Vector etai = Hi*rji;
            const Scalar Wi = kernel(etai.magnitude(), 1.0)/W0;
            CHECK(Wi >= 0.0 and Wi <= 1.0);
            weightSum += Wi;
            phii += Wi*phiij/(phiij*phiij + tiny) * SymTensor::one();
            phimin = min(phimin, Wi*phiij + (1.0 - Wi)*phimin);
          });

          // Increment the overall limiter.
          phi = (phii*phi).Symmetric();
        }

        // Final safety.
        phi *= phimin;
        result(fieldi, i) = phi;
    }
  }

  return result;
}

}
