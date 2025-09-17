//---------------------------------Spheral++----------------------------------//
// CylinderSolidBoundary -- cylinder with finite length solid boundary for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/CylinderSolidBoundary.hh"

#include <string>
using std::string;

namespace Spheral {

template<typename Dimension>
CylinderSolidBoundary<Dimension>::
CylinderSolidBoundary(const Vector& point, 
                      const Vector& axis, 
                      const Scalar radius, 
                      const Scalar length,
                      const RotationType& angularVelocity):
  SolidBoundaryBase<Dimension>(),
  mPoint(point),
  mAxis(axis),
  mRadius(radius),
  mLength(length),
  mVelocity(Vector::zero),
  mAngularVelocity(angularVelocity){
}

template<typename Dimension>
CylinderSolidBoundary<Dimension>::
~CylinderSolidBoundary(){
}

template<typename Dimension>
typename Dimension::Vector
CylinderSolidBoundary<Dimension>::
distance(const Vector& position) const { 
  const auto p = position-mPoint;
  const auto pnMag = p.dot(mAxis);
  const auto pn = pnMag * mAxis;
  const auto paxis = (pnMag > 0 ? max(pnMag-mLength,0.0) : pnMag)*mAxis;
  const auto pr = p-pn;
  return (pr.magnitude() - mRadius)*pr.unitVector() + paxis;
}

template<typename Dimension>
typename Dimension::Vector
CylinderSolidBoundary<Dimension>::
localVelocity(const Vector& position) const { 
  // Calculate the vector from the axis of rotation to the position
  const auto p = position - mPoint;
  const auto pnMag = p.dot(mAxis);
  const auto pn = pnMag * mAxis;
  const auto r = p - pn; // Radial vector from the axis to the position

  // Calculate the tangential velocity due to angular velocity
  const auto tangentialVelocity = DEMDimension<Dimension>::cross(r,mAngularVelocity);
  return mVelocity + tangentialVelocity;
}

template<typename Dimension>
void
CylinderSolidBoundary<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {   
  const auto boundaryKey = "CylinderSolidBoundary_" + std::to_string(std::abs(this->uniqueIndex()));
  const auto pointKey = boundaryKey +"_point";
  const auto axisKey = boundaryKey +"_axis";
  const auto radiusKey = boundaryKey +"_radius";
  const auto lengthKey = boundaryKey +"_length";
  const auto velocityKey = boundaryKey +"_velocity";
  const auto angularVelocityKey = boundaryKey + "_angularVelocity"; // Key for angular velocity
  state.enroll(pointKey,mPoint);
  state.enroll(axisKey,mAxis);
  state.enroll(radiusKey,mRadius);
  state.enroll(lengthKey,mLength);
  state.enroll(velocityKey,mVelocity);
  state.enroll(angularVelocityKey, mAngularVelocity); // Enroll angular velocity
}

template<typename Dimension>
void
CylinderSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mPoint += multiplier*mVelocity;
  // If we want more complex rotation we'll need more complex logic here.
}

//------------------------------------------------------------------------------
// Restart
//------------------------------------------------------------------------------
template<typename Dimension>
void
CylinderSolidBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mPoint, pathName + "/point");
  file.write(mAxis, pathName + "/axis");
  file.write(mRadius, pathName + "/radius");
  file.write(mLength, pathName + "/length");
  file.write(mVelocity, pathName + "/velocity");
  file.write(mAngularVelocity, pathName + "/omega"); // Write angular velocity
}


template<typename Dimension>
void
CylinderSolidBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mPoint, pathName + "/point");
  file.read(mAxis, pathName + "/axis");
  file.read(mRadius, pathName + "/radius");
  file.read(mLength, pathName + "/length");
  file.read(mVelocity, pathName + "/velocity");
  file.read(mAngularVelocity, pathName + "/omega"); // Read angular velocity
}


}