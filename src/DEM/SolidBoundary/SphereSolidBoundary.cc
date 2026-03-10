//---------------------------------Spheral++----------------------------------//
// SphereSolidBoundary -- N-dimensional spherical solid boundary for DEM.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/SphereSolidBoundary.hh"

#include <cmath>
#include <string>
using std::string;

namespace Spheral {

template<typename Dimension>
SphereSolidBoundary<Dimension>::
SphereSolidBoundary(const Vector& center,
                    const Scalar  radius,
                    const RotationType& angularVelocity):
  SolidBoundaryBase<Dimension>(),
  mCenter(center),
  mRadius(radius),
  mVelocity(Vector::zero()),
  mAngularVelocity(angularVelocity){
}

template<typename Dimension>
SphereSolidBoundary<Dimension>::
~SphereSolidBoundary(){
}


template<typename Dimension>
typename Dimension::Vector
SphereSolidBoundary<Dimension>::
distance(const Vector& position) const { 
  const auto p = position - mCenter;
  return p - p.unitVector()*mRadius;
}

template<typename Dimension>
typename Dimension::Vector
SphereSolidBoundary<Dimension>::
localVelocity(const Vector& position) const { 
  const auto rVector = (position - mCenter).unitVector()*mRadius;
  return mVelocity + DEMDimension<Dimension>::cross(mAngularVelocity,rVector);
}

template<typename Dimension>
void
SphereSolidBoundary<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  const auto boundaryKey = "SphereSolidBoundary_" + std::to_string(std::abs(this->uniqueIndex()));
  const auto pointKey = boundaryKey +"_point";
  const auto velocityKey = boundaryKey +"_velocity";

  state.enroll(pointKey,mCenter);
  state.enroll(pointKey,mVelocity);

}

template<typename Dimension>
void
SphereSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mCenter += multiplier*mVelocity;
}

//------------------------------------------------------------------------------
// Restart
//------------------------------------------------------------------------------
template<typename Dimension>
void
SphereSolidBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mAngularVelocity, pathName + "/omega");
  file.write(mCenter, pathName + "/center");
  file.write(mRadius, pathName + "/radius");
  file.write(mVelocity, pathName + "/velocity");
}


template<typename Dimension>
void
SphereSolidBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mAngularVelocity, pathName + "/omega");
  file.read(mCenter, pathName + "/center");
  file.read(mRadius, pathName + "/radius");
  file.read(mVelocity, pathName + "/velocity");
}

}