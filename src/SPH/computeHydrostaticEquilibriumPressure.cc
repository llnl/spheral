//------------------------------------------------------------------------------
// Compute the SPH pressure per point that will produce the input accelerations
// given our standard SPH momentum equation.
// We only cover the 3D case here.
//------------------------------------------------------------------------------
#include "Eigen/Core"
#include "Eigen/Sparse"
#include "Eigen/SVD"
#include "Eigen/LU"

#include "Geometry/Dimension.hh"
#include "computeHydrostaticEquilibriumPressure.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/globalNodeIDs.hh"

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

void
computeSPHHydrostaticEquilibriumPressure(const DataBase<Dim<3> >& db,
                                         const TableKernel<Dim<3> >& W,
                                         const FieldList<Dim<3>, Dim<3>::Vector>& acceleration,
                                         FieldList<Dim<3>, Dim<3>::Scalar>& pressure) {

  // Pre-conditions.
  const unsigned numNodeLists = db.numFluidNodeLists();
  REQUIRE(acceleration.size() == numNodeLists);
  REQUIRE(pressure.size() == numNodeLists);

  typedef Dim<3> Dimension;
  typedef Dim<3>::Scalar Scalar;
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::SymTensor SymTensor;

  // Grab the state fields from the data base.
  const FieldList<Dimension, Vector> position = db.fluidPosition();
  const FieldList<Dimension, Scalar> mass = db.fluidMass();
  const FieldList<Dimension, Scalar> massDensity = db.fluidMassDensity();
  const FieldList<Dimension, SymTensor> H = db.fluidHfield();
  CHECK(position.size() == numNodeLists);
  CHECK(mass.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(position.numGhostElements() == 0);

  // Zero out the result.
  pressure = 0.0;

  // We assume that the Connectivity in the DataBase has already been updated!
  const ConnectivityMap<Dim<3> >& connectivityMap = db.connectivityMap();

  // Get the global IDs for all nodes.
  const FieldList<Dim<3>, size_t> globalIDs = globalNodeIDs<Dim<3>, DataBase<Dim<3> >::ConstFluidNodeListIterator>(db.fluidNodeListBegin(), db.fluidNodeListEnd());

  // Build the sparse matrix that represents the full pressure gradient operator.
  // We have one of these matrix operators for each dimension, hence the 3 vector.
  const unsigned n = pressure.numInternalElements();
  vector<Eigen::SparseMatrix<double, Eigen::RowMajor> > M(3, Eigen::SparseMatrix<double, Eigen::RowMajor>(n, n));
  vector<Eigen::VectorXd> s(3, Eigen::VectorXd(n));
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = pressure[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const auto iglobal = globalIDs(nodeListi, i);

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar& rhoi = massDensity(nodeListi, i);
      const Vector& gi = acceleration(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();

      // Build the s vector.
      s[0][iglobal] = rhoi*gi.x();
      s[1][iglobal] = rhoi*gi.y();
      s[2][iglobal] = rhoi*gi.z();

      // Get the neighbors for this node (in this NodeList).
      const vector<vector<int> >& connectivity = connectivityMap.connectivityForNode(nodeListi, i);
      CHECK(connectivity.size() == numNodeLists);
      for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        for (vector<int>::const_iterator jItr = connectivity[nodeListj].begin();
             jItr != connectivity[nodeListj].end();
             ++jItr) {
          const unsigned j = *jItr;
          const auto jglobal = globalIDs(nodeListj, j);

          // State for node j.
          const Vector& rj = position(nodeListj, j);
          const Scalar& mj = mass(nodeListj, j);
          const Scalar& rhoj = massDensity(nodeListj, j);
          const SymTensor& Hj = H(nodeListj, j);
          const Scalar Hdetj = Hj.Determinant();

          // Kernel gradient.
          const Vector rji = rj - ri;
          const Vector etai = Hi*rji;
          const Vector etaj = Hj*rji;
          const Scalar etaMagi = etai.magnitude();
          const Scalar etaMagj = etaj.magnitude();
          const Vector gradWij = 0.5*(etai*W.gradValue(etaMagi, Hdeti) +
                                      etaj*W.gradValue(etaMagj, Hdetj));
          M[0].insert(iglobal, jglobal) = mj/rhoj*gradWij.x();
          M[1].insert(iglobal, jglobal) = mj/rhoj*gradWij.y();
          M[2].insert(iglobal, jglobal) = mj/rhoj*gradWij.z();
        }
      }
    }
  }

  // 
}

}
