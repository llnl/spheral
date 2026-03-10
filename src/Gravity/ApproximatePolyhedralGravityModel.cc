//---------------------------------Spheral++----------------------------------//
// Approximate Polhedral Model
///----------------------------------------------------------------------------//
//   Model used to calculate the gravitational fields of an irregular shaped 
//   body by apply numerical quadrature to a discetization of its surface.
//
//   Pearl, Hitt, A Fast Quadrature-Based Gravity Model for the Homogeneous 
//   Polyhedron, MNRAS, 492(1), 2020.
//----------------------------------------------------------------------------//
// J. M. Pearl 2022
//----------------------------------------------------------------------------//

#include "ApproximatePolyhedralGravityModel.hh"

namespace Spheral {

//--------------------------------------------------------------------------------------------
// Constructor 
//--------------------------------------------------------------------------------------------
ApproximatePolyhedralGravityModel::
ApproximatePolyhedralGravityModel(const GeomPolyhedron & poly, const double Mass, const double G):
  mNumQuadraturePoints(poly.facets().size()),
  mQuadraturePoints(poly.facetCentroids()),
  mValues(poly.facetAreaVectors()),
  mResolutions(poly.facets().size(), std::numeric_limits<double>::max()) {
    
  auto rho = Mass/poly.volume();
  const auto& facets = poly.facets();
  const auto& verts = poly.vertices();
  for (auto i = 0u; i < facets.size(); ++i) {
    const auto centroid = facets[i].position();
    const auto& iverts = facets[i].ipoints();
    const auto nverts = iverts.size();
    for (auto j = 0u; j < nverts; ++j) {
      mResolutions[i] = std::min(mResolutions[i], (0.5*(verts[iverts[j]] + verts[iverts[(j+1u)%nverts]]) - centroid).magnitude());
    }
    mValues[i] *= G*rho;
  }
}

//--------------------------------------------------------------------------------------------
// Destructor 
//--------------------------------------------------------------------------------------------
ApproximatePolyhedralGravityModel::
~ApproximatePolyhedralGravityModel() {
}

//--------------------------------------------------------------------------------------------
// acceleration
//--------------------------------------------------------------------------------------------
Dim<3>::Vector
ApproximatePolyhedralGravityModel::
acceleration(const Dim<3>::Vector& position) const {
  
  Dim<3>::Vector acceleration = Vector::zero();

  const std::vector<Vector>& GrhoAn = this-> values();
  const std::vector<Vector>& quadPoints = this->quadraturePoints();
  const std::vector<Scalar>& res = this->resolutions();

  for(unsigned int i=0; i < this->numQuadraturePoints(); i++){
    Scalar r = (quadPoints[i] - position).magnitude();
    acceleration -= GrhoAn[i]/std::max(r,res[i]);
  }

  return acceleration;
}

//--------------------------------------------------------------------------------------------
// potential
//--------------------------------------------------------------------------------------------
Dim<3>::Scalar
ApproximatePolyhedralGravityModel::
potential(const Dim<3>::Vector& position) const {
  
  Scalar potential = 0.0;

  const std::vector<Vector>& GrhoAn = this-> values();
  const std::vector<Vector>& quadPoints = this->quadraturePoints();
  const std::vector<Scalar>& res = this->resolutions();

  for(unsigned int i=0; i < this->numQuadraturePoints(); i++){
    Vector r = (quadPoints[i] - position);
    potential -= GrhoAn[i].dot(r)/std::max(r.magnitude(),0.1*res[i]);
  }

  return 0.5*potential;
}

}


