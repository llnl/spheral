//---------------------------------Spheral++----------------------------------//
// A finite-volume based viscosity.  Assumes you have constructred the 
// tessellation in the state.
//
// Created by JMO, Tue Aug 13 09:43:37 PDT 2013
//----------------------------------------------------------------------------//
#include "FiniteVolumeViscosity.hh"
#include "DataOutput/Restart.hh"
#include "Field/FieldList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "NodeList/FluidNodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Material/EquationOfState.hh"
#include "Boundary/Boundary.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "Mesh/Mesh.hh"
#include "Utilities/Timer.hh"

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
// Construct with the given value for the linear and quadratic coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
FiniteVolumeViscosity<Dimension>::
FiniteVolumeViscosity(const Scalar Clinear,
                      const Scalar Cquadratic,
                      const TableKernel<Dimension>& WT):
  ArtificialViscosityView<Dimension, Scalar>(Clinear, Cquadratic, WT) {
}

//------------------------------------------------------------------------------
// Main method -- compute the QPi (P/rho^2) artificial viscosity
//------------------------------------------------------------------------------
template<typename Dimension>
void
FiniteVolumeViscosity<Dimension>::
QPiij(Scalar& QPiij, Scalar& QPiji,      // result for QPi (Q/rho^2)
      Scalar& Qij, Scalar& Qji,          // result for viscous pressure
      const unsigned nodeListi, const unsigned i, 
      const unsigned nodeListj, const unsigned j,
      const Vector& xi,
      const SymTensor& Hi,
      const Vector& etai,
      const Vector& vi,
      const Scalar rhoi,
      const Scalar csi,
      const Vector& xj,
      const SymTensor& Hj,
      const Vector& etaj,
      const Vector& vj,
      const Scalar rhoj,
      const Scalar csj,
      const FieldList<Dimension, Scalar>& fCl,
      const FieldList<Dimension, Scalar>& fCq,
      const FieldList<Dimension, Tensor>& DvDx) const {

  // Preconditions
  REQUIRE(fCl.size() == fCq.size());
  REQUIRE((not mBalsaraShearCorrection) or DvDx.size() > std::max(nodeListi, nodeListj));

  // A few useful constants
  const auto multipliers = fCl.size() > 0u;

  // Find our linear and quadratic coefficients
  const auto fCli = multipliers ? fCl(nodeListi, i) : 1.0;
  const auto fCqi = multipliers ? fCq(nodeListi, i) : 1.0;
  const auto fClj = multipliers ? fCl(nodeListj, j) : 1.0;
  const auto fCqj = multipliers ? fCq(nodeListj, j) : 1.0;
  const auto fshear = (mBalsaraShearCorrection ?
                       0.5*(this->calcBalsaraShearCorrection(DvDx(nodeListi, i), Hi, csi) +
                            this->calcBalsaraShearCorrection(DvDx(nodeListj, j), Hj, csj)) :
                       1.0);
  const auto Clij = 0.5*(fCli + fClj)*fshear * mClinear;
  const auto Cqij = 0.5*(fCqi + fCqj)*fshear * mCquadratic;

  // Compute the pair QPi
  const auto xji = xj - xi;
  const auto xjihat = xji.unitVector();
  const auto hi = 1.0/(Hi*xjihat).magnitude();
  const auto hj = 1.0/(Hj*xjihat).magnitude();
  const auto DvDxi = min(0.0, DvDx(nodeListi, i).Trace());
  const auto DvDxj = min(0.0, DvDx(nodeListj, j).Trace());
  QPiij = (-Clij*csi*DvDxi + Cqij*hi*DvDxi*DvDxi)*hi/rhoi;
  QPiji = (-Clij*csj*DvDxj + Cqij*hj*DvDxj*DvDxj)*hj/rhoj;
  Qij = rhoi*rhoi*QPiij;
  Qji = rhoi*rhoi*QPiji;
}

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
      const vector<int>& faces = zonei.faceIDs();
      Vi = zonei.volume();
      for (vector<int>::const_iterator fitr = faces.begin();
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
