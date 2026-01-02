//---------------------------------Spheral++----------------------------------//
// FieldListFunctions -- A set of global functions which can be applied to
// FieldLists.
//
// Created by JMO, Wed Dec  6 21:09:29 PST 2000
//----------------------------------------------------------------------------//
#include "FieldListFunctions.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
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
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, GradientType>("grad", (*fieldItr)->nodeList()));
    flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Loop over all the elements in the input FieldList.
  for (InternalNodeIterator<Dimension> nodeItr = fieldList.internalNodeBegin();
       nodeItr < fieldList.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // We will do the batch of master nodes associated with this node together.
      // Set the neighbor information.
      vector<vector<int>> masterLists, coarseNeighbors, refineNeighbors;
      fieldList.setMasterNodeLists(position(nodeItr), Hfield(nodeItr), masterLists, coarseNeighbors);

      // Now loop over all the master nodes.
      for (MasterNodeIterator<Dimension> masterItr = fieldList.masterNodeBegin(masterLists);
           masterItr < fieldList.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);

        // Set the refined neighbor information for this master node.
        fieldList.setRefineNodeLists(position(masterItr), Hfield(masterItr), coarseNeighbors, refineNeighbors);

        // Loop over the refined neighbors.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& rhoi = rho(masterItr);
        const DataType& fieldi = fieldList(masterItr);

        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin(refineNeighbors);
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {

          const Vector& rj = position(neighborItr);
          const SymTensor& Hj = Hfield(neighborItr);
          const Scalar& mj = mass(neighborItr);
          const DataType& fieldj = fieldList(neighborItr);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Vector etaiNorm = etai.unitVector();
          const Vector etajNorm = etaj.unitVector();

          // Get the symmetrized kernel gradient for this node pair.
          Vector gradWij;
          switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
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
          result(masterItr) += mj*(fieldj - fieldi)*gradWij;
        }

        // Normalize by density.
        CHECK(rhoi > 0.0);
        result(masterItr) /= rhoi;

        // This master node is finished.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (typename vector< vector<bool> >::const_iterator
         flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
    if (checkcount > 0) {
      cerr << "Error in FieldList::gradient: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
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
  const auto vectorSize = (fieldList.numInternalElements() > 0 ?
                           fieldList(fieldList.internalNodeBegin()).size() :
                           0);
  
  // Return FieldList.
  FieldList<Dimension, std::vector<GradientType>> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, std::vector<DataType>>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, std::vector<GradientType>>("grad",
                                                                   (*fieldItr)->nodeList(),
                                                                   std::vector<GradientType>(vectorSize, GradientType::zero())));
    flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Loop over all the elements in the input FieldList.
  for (InternalNodeIterator<Dimension> nodeItr = fieldList.internalNodeBegin();
       nodeItr < fieldList.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // We will do the batch of master nodes associated with this node together.
      // Set the neighbor information.
      vector<vector<int>> masterLists, coarseNeighbors, refineNeighbors;
      fieldList.setMasterNodeLists(position(nodeItr), Hfield(nodeItr), masterLists, coarseNeighbors);

      // Now loop over all the master nodes.
      for (MasterNodeIterator<Dimension> masterItr = fieldList.masterNodeBegin(masterLists);
           masterItr < fieldList.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);

        // Set the refined neighbor information for this master node.
        fieldList.setRefineNodeLists(position(masterItr), Hfield(masterItr), coarseNeighbors, refineNeighbors);

        // Loop over the refined neighbors.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& rhoi = rho(masterItr);
        const std::vector<DataType>& fieldi = fieldList(masterItr);

        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin(refineNeighbors);
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {

          const Vector& rj = position(neighborItr);
          const SymTensor& Hj = Hfield(neighborItr);
          const Scalar& mj = mass(neighborItr);
          const std::vector<DataType>& fieldj = fieldList(neighborItr);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Vector etaiNorm = etai.unitVector();
          const Vector etajNorm = etaj.unitVector();

          // Get the symmetrized kernel gradient for this node pair.
          Vector gradWij;
          switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
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
            CHECK(result(masterItr).size() == vectorSize &&
                  fieldj.size() == vectorSize &&
                  fieldi.size() == vectorSize);
            result(masterItr)[m] += mj*(fieldj[m] - fieldi[m])*gradWij;
          }
        }

        // Normalize by density.
        CHECK(rhoi > 0.0);
        for (auto m = 0u; m < vectorSize; ++m) {
          CHECK(result(masterItr).size() == vectorSize);
          result(masterItr)[m] /= rhoi;
        }

        // This master node is finished.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (typename vector< vector<bool> >::const_iterator
         flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
    if (checkcount > 0) {
      cerr << "Error in FieldList::gradient: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
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
  result = GradientType::zero;
  
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
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator
         fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, SymTensor>("limiter", (*fieldItr)->nodeList()));
    flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Loop over all the elements in the input FieldList.
  for (InternalNodeIterator<Dimension> nodeItr = fieldList.internalNodeBegin();
       nodeItr != fieldList.internalNodeEnd();
       ++nodeItr) {

    // Check if this node has been done yet.
    if (!flagNodeDone[nodeItr.fieldID()][nodeItr.nodeID()]) {

      // We will do the batch of master nodes associated with this node together.
      // Set the neighbor information.
      vector<vector<int>> masterLists, coarseNeighbors, refineNeighbors;
      fieldList.setMasterNodeLists(position(nodeItr), Hfield(nodeItr), masterLists, coarseNeighbors);

      // Now loop over all the master nodes.
      for (MasterNodeIterator<Dimension> masterItr = fieldList.masterNodeBegin(masterLists);
           masterItr < fieldList.masterNodeEnd();
           ++masterItr) {
        CHECK(flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] == false);

        // Set the refined neighbor information for this master node.
        fieldList.setRefineNodeLists(position(masterItr), Hfield(masterItr), coarseNeighbors, refineNeighbors);

        // State for node i.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const DataType& fieldi = fieldList(masterItr);
        const GradientType& gradi = gradient(masterItr);

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
          for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin(refineNeighbors);
               neighborItr != fieldList.refineNodeEnd();
               ++neighborItr) {

            // State for node j.
            const Vector& rj = position(neighborItr);
            //const SymTensor& Hj = Hfield(neighborItr);
            const DataType& fieldj = fieldList(neighborItr);

            // Compute the pair-wise limiting needed.
            const Vector rji = rj - ri;
            const DataType dF = fieldj - fieldi;
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
          }

          // Increment the overall limiter.
          phi = (phii*phi).Symmetric();
        }

        // Final safety.
        phi *= phimin;
        result(masterItr) = phi;
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  return result;
}

}
