//---------------------------------Spheral++----------------------------------//
// FieldListFunctionsMash -- A set of global functions which can be applied to
// FieldLists using MASH prescriptions.
//
// Created by JMO, Wed Dec  6 21:09:29 PST 2000
//----------------------------------------------------------------------------//
#include "FieldListFunctionsMash.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
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
gradientMash(const FieldList<Dimension, DataType>& fieldList,
             const FieldList<Dimension, typename Dimension::Vector>& position,
             const FieldList<Dimension, typename Dimension::Scalar>& weight,
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
  for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, GradientType>("grad mash", (*fieldItr)->nodeList()));
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
        const Scalar& weighti = weight(masterItr);
        const DataType& Fi = fieldList(masterItr);

        Scalar normalization = weighti*kernel(0.0, 1.0);
        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin(refineNeighbors);
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {
          if (neighborItr != masterItr) {

            const Vector& rj = position(neighborItr);
            const SymTensor& Hj = Hfield(neighborItr);
            const Scalar& weightj = weight(neighborItr);
            const DataType& Fj = fieldList(neighborItr);

            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;
            const Scalar etaiMag = etai.magnitude();
            const Scalar etajMag = etaj.magnitude();
            const Vector etaiNorm = etai.unitVector();
            const Vector etajNorm = etaj.unitVector();

            // Get the symmetrized kernel gradient for this node pair.
            Scalar Wij;
            Vector gradWij;
            switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
            case NeighborSearchType::GatherScatter:
              Wij = 0.5*(kernel(etaiMag, 1.0) + kernel(etajMag, 1.0));
              gradWij = 0.5*(Hi*etaiNorm*kernel.grad(etaiMag, 1.0) + 
                             Hj*etajNorm*kernel.grad(etajMag, 1.0));
              break;

            case NeighborSearchType::Gather:
              Wij = kernel(etaiMag, 1.0);
              gradWij = Hi*etaiNorm*kernel.grad(etaiMag, 1.0);
              break;

            case NeighborSearchType::Scatter:
              Wij = kernel(etajMag, 1.0);
              gradWij = Hj*etajNorm*kernel.grad(etajMag, 1.0);
              break;

            default:
              VERIFY2(false, "Unhandled neighbor search type.");
            }

            // Add this nodes contribution to the master value.
            normalization += weightj*Wij;
            result(masterItr) += weightj*(Fj - Fi)*gradWij;
          }
        }

        CHECK(normalization > 0.0);
        result(masterItr) /= normalization;

        // This master node is finished.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (typename vector< vector<bool> >::const_iterator flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
    if (checkcount > 0) {
      cerr << "Error in FieldList::smoothFieldsMash: Not all values determined on exit "
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
FieldList<Dimension, typename MathTraits<Dimension, DataType>::GradientType>
gradientMash2(const FieldList<Dimension, DataType>& fieldList,
              const FieldList<Dimension, typename Dimension::Vector>& position,
              const FieldList<Dimension, typename Dimension::Scalar>& weight,
              const FieldList<Dimension, typename Dimension::Scalar>& weightDensity,
              const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
              const TableKernel<Dimension>& kernel) {

  // Typedef's to ease typing/understandability.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename MathTraits<Dimension, DataType>::GradientType GradientType;

  // Return FieldList.
  FieldList<Dimension, GradientType> result;
  FieldList<Dimension, DataType> smoothFieldList;
  FieldList<Dimension, Scalar> a;
  FieldList<Dimension, Scalar> b;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  smoothFieldList.copyFields();
  a.copyFields();
  b.copyFields();
  for (typename FieldList<Dimension, DataType>::const_iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, GradientType>("gradient", (*fieldItr)->nodeList()));
    smoothFieldList.appendField(Field<Dimension, DataType>("smooth", (*fieldItr)->nodeList()));
    a.appendField(Field<Dimension, Scalar>("a", (*fieldItr)->nodeList()));
    b.appendField(Field<Dimension, Scalar>("b", (*fieldItr)->nodeList()));
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
        //const Scalar& weighti = weight(masterItr);
        const Scalar& densityi = weightDensity(masterItr);
        //const DataType& fieldi = fieldList(masterItr);

        // First identify the fitting parameters for a linear model (w_j = a + b eta_j).
        Scalar totalWeight = 0.0;
        Vector rCenter;
        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin(refineNeighbors);
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {

          const Vector& rj = position(neighborItr);
          const SymTensor& Hj = Hfield(neighborItr);
          const Scalar& weightj = weight(neighborItr);
          //const Scalar& densityj = weightDensity(neighborItr);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Scalar etaiMag = etai.magnitude();
          const Scalar etajMag = etaj.magnitude();

          // Get the symmetrized kernel weighting for this node pair.
          Scalar Wij;
          switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
          case NeighborSearchType::GatherScatter:
            Wij = 0.5*(kernel(etaiMag, 1.0) + 
                       kernel(etajMag, 1.0));
            break;

          case NeighborSearchType::Gather:
            Wij = kernel(etaiMag, 1.0);
            break;

          case NeighborSearchType::Scatter:
            Wij = kernel(etajMag, 1.0);
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Sum this nodes contributions to the fitting parameters.
          totalWeight += weightj*Wij;
          rCenter -= weightj*Wij*rij;
        }

        // Finally we can calculate the correction factors.
        CHECK(totalWeight > 0.0);
        rCenter /= totalWeight;
        a(masterItr) = densityi;
        b(masterItr) = 2.0*densityi*rCenter.x()*Hi.xx()*Hi.xx();
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (typename vector< vector<bool> >::const_iterator flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
    if (checkcount > 0) {
      cerr << "Error in FieldList::gradientMash2: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  for (typename vector< vector<bool> >::iterator flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    for (typename vector<bool>::iterator itr = flagNodeItr->begin();
         itr < flagNodeItr->end();
         ++itr) {
      *itr = false;
    }
  }

  // Begin by finding a smoothed estimate of the FieldList.
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
        //const Scalar& weighti = weight(masterItr);

        // First calculate the smoothed estimate of Fi.
        Scalar normalization = 0.0;
        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin(refineNeighbors);
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {

          const Vector& rj = position(neighborItr);
          const SymTensor& Hj = Hfield(neighborItr);
          const Scalar& weightj = weight(neighborItr);
          const DataType& Fj = fieldList(neighborItr);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Scalar etaiMag = etai.magnitude();
          const Scalar etajMag = etaj.magnitude();

          // Get the symmetrized kernel gradient for this node pair.
          Scalar Wij;
          switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
          case NeighborSearchType::GatherScatter:
            Wij = 0.5*(kernel(etaiMag, 1.0) +
                       kernel(etajMag, 1.0))*(a(masterItr) + b(masterItr)*rij.x());
            break;

          case NeighborSearchType::Gather:
            Wij = kernel(etaiMag, 1.0)*(a(masterItr) + b(masterItr)*rij.x());
            break;

          case NeighborSearchType::Scatter:
            Wij = kernel(etajMag, 1.0)*(a(masterItr) + b(masterItr)*rij.x());
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }

          // Add this nodes contribution to the master value.
          normalization += weightj*Wij;
          smoothFieldList(masterItr) += Fj*weightj*Wij;

        }
        CHECK(normalization > 0.0);
        smoothFieldList(masterItr) /= normalization;

        // This master node is finished.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (typename vector< vector<bool> >::const_iterator flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
    if (checkcount > 0) {
      cerr << "Error in FieldList::gradientMash2: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  // Reset the flag node lists.
  for (typename vector< vector<bool> >::iterator outerItr = flagNodeDone.begin();
       outerItr < flagNodeDone.end();
       ++outerItr) {
    for (typename vector<bool>::iterator innerItr = outerItr->begin();
       innerItr < outerItr->end();
       ++innerItr) {
      *innerItr = false;
    }
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
        //const Scalar& weighti = weight(masterItr);
        //const DataType& Fi = fieldList(masterItr);
        //const DataType& Fi = smoothFieldList(masterItr);

        Scalar normalization = 0.0;
        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin(refineNeighbors);
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {

          const Vector& rj = position(neighborItr);
          const SymTensor& Hj = Hfield(neighborItr);
          const Scalar& weightj = weight(neighborItr);
          const DataType& Fj = fieldList(neighborItr);

          const Vector rij = ri - rj;
          const Vector etai = Hi*rij;
          const Vector etaj = Hj*rij;
          const Scalar etaiMag = etai.magnitude();
          const Scalar etajMag = etaj.magnitude();
          const Vector etaiNorm = etai.unitVector();
          const Vector etajNorm = etaj.unitVector();

          // Get the symmetrized kernel gradient for this node pair.
          Scalar Wij;
          Vector gradWij;
          switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
          case NeighborSearchType::GatherScatter:
            Wij = 0.5*(kernel(etaiMag, 1.0) + kernel(etajMag, 1.0));
            gradWij = 0.5*(Hi*etaiNorm*kernel.grad(etaiMag, 1.0) + 
                           Hj*etajNorm*kernel.grad(etajMag, 1.0));
            break;

          case NeighborSearchType::Gather:
            Wij = kernel(etaiMag, 1.0);
            gradWij = Hi*etaiNorm*kernel.grad(etaiMag, 1.0);
            break;

          case NeighborSearchType::Scatter:
            Wij = kernel(etajMag, 1.0);
            gradWij = Hj*etajNorm*kernel.grad(etajMag, 1.0);
            break;

          default:
            VERIFY2(false, "Unhandled neighbor search type.");
          }
          Wij *= a(masterItr) + b(masterItr)*rij.x();
          gradWij = ((a(masterItr) + b(masterItr)*rij.x())*gradWij -
                     Vector(b(masterItr)*Wij));

          // Add this nodes contribution to the master value.
          normalization += weightj*Wij;
          result(masterItr) += weightj*(Fj - smoothFieldList(masterItr))*gradWij;
        }

        CHECK(normalization > 0.0);
        result(masterItr) /= normalization;

        // This master node is finished.
        flagNodeDone[masterItr.fieldID()][masterItr.nodeID()] = true;
      }
    }
  }

  // After we're done, all nodes in all NodeLists should be flagged as done.
  for (typename vector< vector<bool> >::const_iterator flagNodeItr = flagNodeDone.begin();
       flagNodeItr < flagNodeDone.end();
       ++flagNodeItr) {
    int checkcount = count(flagNodeItr->begin(), flagNodeItr->end(), false);
    if (checkcount > 0) {
      cerr << "Error in FieldList::smoothFieldsMash: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  return result;
}

}
