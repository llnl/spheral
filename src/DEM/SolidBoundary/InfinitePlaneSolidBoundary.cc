//---------------------------------Spheral++----------------------------------//
// InfinitePlaneSolidBoundary -- rigid planar wall contact for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/InfinitePlaneSolidBoundary.hh"

#include <cmath>
#include <string>
using std::string;

namespace Spheral {

template<typename Dimension>
InfinitePlaneSolidBoundary<Dimension>::
InfinitePlaneSolidBoundary(const Vector& point, 
                           const Vector& normal, 
                           const RotationType& angularVelocity):
  SolidBoundaryBase<Dimension>(),
  mPoint(point),
  mNormal(normal),
  mVelocity(Vector::zero),
  mAngularVelocity(angularVelocity){
}

template<typename Dimension>
InfinitePlaneSolidBoundary<Dimension>::
~InfinitePlaneSolidBoundary(){
}

template<typename Dimension>
typename Dimension::Vector
InfinitePlaneSolidBoundary<Dimension>::
distance(const Vector& position) const { 
  return (position - mPoint).dot(mNormal)*mNormal;
}

template<typename Dimension>
typename Dimension::Vector
InfinitePlaneSolidBoundary<Dimension>::
localVelocity(const Vector& position) const { 
  return mVelocity + DEMDimension<Dimension>::cross((position-mPoint),mAngularVelocity);  // Include angular velocity effect
}

template<typename Dimension>
void
InfinitePlaneSolidBoundary<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  const auto boundaryKey = "InfinitePlaneSolidBoundary_" + std::to_string(std::abs(this->uniqueIndex()));
  const auto pointKey = boundaryKey +"_point";
  const auto velocityKey = boundaryKey +"_velocity";
  const auto normalKey = boundaryKey +"_normal";
  const auto angularVelocityKey = boundaryKey +"_angularVelocity";
  state.enroll(pointKey,mPoint);
  state.enroll(velocityKey,mVelocity);
  state.enroll(normalKey,mNormal);
  state.enroll(angularVelocityKey, mAngularVelocity);  // Enroll angular velocity
}

template<typename Dimension>
void
InfinitePlaneSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mPoint += multiplier*mVelocity;
  // Update the normal vector based on angular velocity
  // Using a small angle approximation: newNormal = oldNormal + dt * (angularVelocity x oldNormal)
  mNormal += multiplier * DEMDimension<Dimension>::cross(mNormal,mAngularVelocity);
  mNormal = mNormal.unitVector();  // Normalize to maintain it as a unit vector
}

//------------------------------------------------------------------------------
// Restart
//------------------------------------------------------------------------------
template<typename Dimension>
void
InfinitePlaneSolidBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mPoint, pathName + "/point");
  file.write(mNormal, pathName + "/normal");
  file.write(mVelocity, pathName + "/velocity");
  file.write(mAngularVelocity, pathName + "/omega");  // Write angular velocity
}


template<typename Dimension>
void
InfinitePlaneSolidBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mPoint, pathName + "/point");
  file.read(mNormal, pathName + "/normal");
  file.read(mVelocity, pathName + "/velocity");
  file.read(mAngularVelocity, pathName + "/omega");  // Read angular velocity
}

}