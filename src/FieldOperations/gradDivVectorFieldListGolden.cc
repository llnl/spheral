//---------------------------------Spheral++----------------------------------//
// FieldListSecondDerivatives
// A set of experimental methods for evaluating the second derivative of a 
// FieldList.
//
// Created by JMO, Wed Dec 18 22:46:54 PST 2002
//----------------------------------------------------------------------------//

#include "FieldListSecondDerivatives.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"

#include <vector>
using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Calculate the gradient of the divergence of a Vector FieldList.
// Golden rule method.
//------------------------------------------------------------------------------
template<typename Dimension>
FieldList<Dimension, typename Dimension::Vector>
gradDivVectorFieldListGolden
(const FieldList<Dimension, typename Dimension::Vector>& fieldList,
 const FieldList<Dimension, typename Dimension::Vector>& position,
 const FieldList<Dimension, typename Dimension::Scalar>& /*weight*/,
 const FieldList<Dimension, typename Dimension::Scalar>& mass,
 const FieldList<Dimension, typename Dimension::Scalar>& rho,
 const FieldList<Dimension, typename Dimension::SymTensor>& Hfield,
 const TableKernel<Dimension>& kernel) {

  // Some convenient typedefs.
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  // Return FieldList.
  FieldList<Dimension, Vector> result;
  vector< vector<bool> > flagNodeDone(fieldList.numFields());
  result.copyFields();
  for (typename FieldList<Dimension, Vector>::const_iterator fieldItr = fieldList.begin();
       fieldItr < fieldList.end(); 
       ++fieldItr) {
    result.appendField(Field<Dimension, Vector>("grad div " + (*fieldItr)->name(), (*fieldItr)->nodeList()));
    flagNodeDone[fieldItr - fieldList.begin()].resize((*fieldItr)->nodeListPtr()->numInternalNodes(), false);
  }

  // Remember the square of the kernel extent.
  const Scalar cutoff2 = FastMath::pow2(kernel.kernelExtent());

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

        // State for this node.
        const Vector& ri = position(masterItr);
        const SymTensor& Hi = Hfield(masterItr);
        const Scalar& mi = mass(masterItr);
        const Scalar& rhoi = rho(masterItr);
        //const Scalar& weighti = weight(masterItr);
        const Vector& fieldi = fieldList(masterItr);
        CONTRACT_VAR(mi);
        CHECK(Hi.Determinant() > 0.0);
        CHECK(mi > 0.0);
        CHECK(rhoi > 0.0);

        // Temp variables to accumulate the grad field and grad rho for this
        // node.
        vector<Tensor> grad2elements(Dimension::nDim);
        Vector gradRho;
        Tensor gradF;

        // Loop over the refined neighbors, and calculate the various
        // sums that contributed to grad div F.
        for (RefineNodeIterator<Dimension> neighborItr = fieldList.refineNodeBegin(refineNeighbors);
             neighborItr < fieldList.refineNodeEnd();
             ++neighborItr) {
          if (neighborItr != masterItr) {

            const Vector& rj = position(neighborItr);
            const SymTensor& Hj = Hfield(neighborItr);
            const Scalar& mj = mass(neighborItr);
            const Scalar& rhoj = rho(neighborItr);
            //const Scalar& weightj = weight(neighborItr);
            const Vector& fieldj = fieldList(neighborItr);
            CHECK(Hj.Determinant() > 0.0);
            CHECK(mj > 0.0);
            CHECK(rhoj > 0.0);

            const Vector rij = ri - rj;
            const Vector etai = Hi*rij;
            const Vector etaj = Hj*rij;
            const Scalar etaiMag = etai.magnitude();
            const Scalar etajMag = etai.magnitude();

            if (etai.magnitude2() < cutoff2 ||
                etaj.magnitude2() < cutoff2) {
              const Vector etaiNorm = etai.unitVector();
              const Vector etajNorm = etaj.unitVector();
              const Vector Hetai = Hi*etaiNorm;
              const Vector Hetaj = Hj*etajNorm;
              const Scalar Hdeti = Hi.Determinant();
              const Scalar Hdetj = Hj.Determinant();

              const Scalar getai = kernel.grad(etaiMag, Hdeti);
              const Scalar getaj = kernel.grad(etajMag, Hdetj);
              const Scalar g2etai = kernel.grad2(etaiMag, Hdeti);
              const Scalar g2etaj = kernel.grad2(etajMag, Hdetj);

              const Vector gWi = Hetai*getai;
              const Vector gWj = Hetaj*getaj;

              const Tensor acki = Hetai.dyad(Hetai);
              const Tensor ackj = Hetaj.dyad(Hetaj);

              const Tensor H2i = Hi*Hi;
              const Tensor H2j = Hj*Hj;

              const Tensor g2Wi = acki*g2etai + 
                (H2i.Transpose() - acki)/(etai.magnitude() + 1.0e-30)*getai;
              const Tensor g2Wj = ackj*g2etaj +
                (H2j.Transpose() - ackj)/(etaj.magnitude() + 1.0e-30)*getaj;

              // Get the symmetrized kernel weighting for this node pair.
              Vector gWij;
              Tensor g2Wij;
              switch((*fieldList.begin())->nodeListPtr()->neighbor().neighborSearchType()) {
              case NeighborSearchType::GatherScatter:
                gWij = 0.5*(gWi + gWj);
                g2Wij = 0.5*(g2Wi + g2Wj);
                break;

              case NeighborSearchType::Gather:
                gWij = gWi;
                g2Wij = g2Wi;
                break;

              case NeighborSearchType::Scatter:
                gWij = gWj;
                g2Wij = g2Wj;
                break;

              default:
                VERIFY2(false, "Unhandled neighbor search type.");
              }

              // Sum the contributions.

//               const Scalar weightij = 0.5*(weighti + weightj);
//               const Vector fji = weightij*(fieldj - fieldi);
//               for (int gamma = 0; gamma < Dimension::nDim; ++gamma) {
//                 for (int beta = 0; beta < Dimension::nDim; ++beta) {
//                   for (int alpha = 0; alpha < Dimension::nDim; ++alpha) {
//                     grad2elements[gamma](beta, alpha) += fji(gamma)*g2Wij(beta, alpha);
//                   }
//                   gradF(gamma, beta) += fji(gamma)*gWij(beta);
//                 }
//                 gradRho(gamma) += weightij*(rhoj - rhoi)*gWij(gamma);
//               }

              const Vector fij = fieldi - fieldj;
              for (int gamma = 0; gamma < Dimension::nDim; ++gamma) {
                for (int beta = 0; beta < Dimension::nDim; ++beta) {
                  for (int alpha = 0; alpha < Dimension::nDim; ++alpha) {
                    grad2elements[gamma](beta, alpha) -= mj*fij(gamma)*g2Wij(beta, alpha);
                  }
                  gradF(gamma, beta) -= mj*fij(gamma)*gWij(beta);
                }
                gradRho(gamma) -= mj*(rhoi - rhoj)*gWij(gamma);
              }
            }
          }
        }

        // Put together the complete derivative.
        CHECK(rhoi > 0.0);
        gradRho /= rhoi;
        gradF /= rhoi;
        for (int beta = 0; beta < Dimension::nDim; ++beta) {
          for (int alpha = 0; alpha < Dimension::nDim; ++alpha) {
            result(masterItr)(beta) += grad2elements[alpha](beta, alpha) -
              gradRho(alpha)*gradF(alpha, beta) -
              gradRho(beta)*gradF(alpha, alpha);
          }
        }
        result(masterItr) /= rhoi;

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
      cerr << "Error in FieldList::smoothFields: Not all values determined on exit "
           << checkcount << endl;
    }
    CHECK(checkcount == 0);
  }

  return result;
}

}
