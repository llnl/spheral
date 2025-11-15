//------------------------------------------------------------------------------
// Centroidally relax a node distribution in a boundary.
// Optionally the user can specify a weighting function for the nodes.
//------------------------------------------------------------------------------
#include "relaxNodeDistribution.hh"
#include "Mesh/Mesh.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Distributed/allReduce.hh"

#include <ctime>
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


template<typename Dimension>
void
relaxNodeDistribution(DataBase<Dimension>& dataBase,
                      const typename Dimension::FacetedVolume& boundary,
                      const std::vector<Boundary<Dimension>*>& /*boundaries*/,
                      const TableKernel<Dimension>& /*W*/,
                      const WeightingFunctor<Dimension>& weightingFunctor,
                      const WeightingFunctor<Dimension>& massDensityFunctor,
                      const double targetMass,
                      const int maxIterations,
                      const double tolerance) {

  typedef Mesh<Dimension> MeshType;
  typedef typename MeshType::Zone Zone;
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Grab the state.
  FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  FieldList<Dimension, Scalar> rho = dataBase.fluidMassDensity();
  FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  FieldList<Dimension, Vector> delta = dataBase.newFluidFieldList(Vector::zero(), "delta");
  const Vector& xmin = boundary.xmin();
  const Vector& xmax = boundary.xmax();
  const double stopTol = tolerance*((xmax - xmin).maxAbsElement());
  const unsigned numNodeLists = dataBase.numNodeLists();

  // Iterate until we either hit the maximum number of iterations or hit the convergence
  // tolerance.
  MeshType mesh;
  double maxDelta = 2.0*stopTol;
  int iter = 0;
  while (iter < maxIterations and maxDelta > stopTol) {
    ++iter;
    delta = Vector::zero();
    maxDelta = 0.0;

    // Read out the node positions to a flat list.
    vector<Vector> generators;
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      copy(position[nodeListi]->begin(), 
           position[nodeListi]->begin() + n,
           back_inserter(generators));
    }

    // Generate the tessellation.
    mesh.reconstruct(generators, boundary);

    // Apply centroidal filtering.
    unsigned k = 0;
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        const Zone& zonei = mesh.zone(k);
        const vector<unsigned>& nodes = zonei.nodeIDs();
        ++k;

        // Compute a weighted centroid.
        double Wsum = 0.0;
        Vector centroid;
        const unsigned nnodes = nodes.size();
        for (unsigned e1 = 0; e1 != nnodes; ++e1) {
          const unsigned e2 = (e1 + 1) % nnodes;
          const Vector n1 = mesh.node(nodes[e1]).position();
          const Vector n2 = mesh.node(nodes[e2]).position();
          const Vector n12 = 0.5*(n1 + n2);
          const double Wj = weightingFunctor(n12, boundary)*(n2 - n1).magnitude();
          Wsum += Wj;
          centroid += Wj*n12;
        }
        CHECK(Wsum > 0.0);
        centroid /= Wsum;

        // Compute the delta to move toward the centroid.
        const Vector& ri = position(nodeListi, i);
        delta(nodeListi, i) = 0.5*(ri + centroid) - ri;

        // If we moved outside the boundary, go back!
        if (not boundary.contains(ri + delta(nodeListi, i))) {
          const Vector p = boundary.closestPoint(ri);
          delta(nodeListi, i) = 0.5*(ri + p) - ri;
        }
        maxDelta = std::max(maxDelta, delta(nodeListi, i).magnitude());
      }
    }
    if (Process::getRank() == 0) cerr << "relaxNodeDistribution iteration " << iter << " maxDelta=" << maxDelta << endl;

    // Apply the change.
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        position(nodeListi, i) += delta(nodeListi, i);
        CHECK(boundary.contains(position(nodeListi, i)));
      }
    }
  }

  // Update the mass and mass density.  We use the last iteration of the mesh 
  // here -- hopefully not too far off.
  unsigned k = 0;
  double Msum = 0.0;
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = position[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Zone& zonei = mesh.zone(k);
      const Vector& ri = position(nodeListi, i);
      const double voli = zonei.volume();
      rho(nodeListi, i) = massDensityFunctor(ri, boundary);
      mass(nodeListi, i) = voli*rho(nodeListi, i);
      Msum += mass(nodeListi, i);
      ++k;
    }
  }
  Msum = allReduce(Msum, SPHERAL_OP_SUM);

  // If needed, rescale the masses.
  if (targetMass > 0.0) {
    CHECK(Msum > 0.0);
    const double f = targetMass/Msum;
    rho *= f;
    mass *= f;
  }
}

}
