//---------------------------------Spheral++----------------------------------//
// A finite-volume based viscosity.  Assumes you have constructred the 
// tessellation in the state.
//
// Created by JMO, Tue Aug 13 09:43:37 PDT 2013
//----------------------------------------------------------------------------//
#include "FiniteVolumeViscosity.hh"
#include "DataBase/State.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Mesh/Mesh.hh"
#include "Utilities/Timer.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Override the base method of computing the velocity gradient with a
// finite-volume approach
//------------------------------------------------------------------------------
template<typename Dimension>
void
FiniteVolumeViscosity<Dimension>::
updateVelocityGradient(const DataBase<Dimension>& dataBase,
                       const State<Dimension>& state,
                       const StateDerivatives<Dimension>& derivs) {
  TIME_BEGIN("FiniteVolumeViscosity_updateVelocityGradient");

  using Zone = typename Mesh<Dimension>::Zone;
  using Face = typename Mesh<Dimension>::Face;

  // Grab the DvDx for updating
  auto DvDx = state.fields(HydroFieldNames::ArtificialViscosityVelocityGradient, Tensor::zero);
  DvDx.Zero();

  // Make a finite-volume estimate of the local (to each Voronoi cell) velocity
  // gradient.
  unsigned nodeListj, j;
  Scalar Vi;
  const Mesh<Dimension>& mesh = state.mesh();
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const unsigned numNodeLists = velocity.numFields();
  for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
    const unsigned n = velocity[nodeListi]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      const Vector& vi = velocity(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      const Zone& zonei = mesh.zone(nodeListi, i);
      const std::vector<int>& faces = zonei.faceIDs();
      Vi = zonei.volume();
      for (std::vector<int>::const_iterator fitr = faces.begin();
           fitr != faces.end();
           ++fitr) {
        const Face& faceij = mesh.face(*fitr);
        const int oppZoneID = faceij.oppositeZoneID(zonei.ID());
        if (Mesh<Dimension>::positiveID(oppZoneID) == Mesh<Dimension>::UNSETID) {
          nodeListj = nodeListi;
          j = i;
        } else {
          mesh.lookupNodeListID(Mesh<Dimension>::positiveID(oppZoneID), nodeListj, j);
        }
        const Vector& vj = velocity(nodeListj, j);
        const Vector vij = 0.5*(vi + vj);
        const Vector dA = faceij.area() * faceij.unitNormal() * sgn(*fitr);
        DvDxi -= vij*dA;
      }
      DvDxi /= Vi;
    }
  }
  TIME_END("FiniteVolumeViscosity_updateVelocityGradient");
}

}
