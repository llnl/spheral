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
               const Scalar length):
  SolidBoundaryBase<Dimension>(),
  mPoint(point),
  mAxis(axis),
  mRadius(radius),
  mLength(length),
  mVelocity(Vector::zero()){
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
  return mVelocity;
}

template<typename Dimension>
void
CylinderSolidBoundary<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {   
  const auto boundaryKey = "CylinderSolidBoundary_" + std::to_string(std::abs(this->uniqueIndex()));
  const auto pointKey = boundaryKey +"_point";
  const auto velocityKey = boundaryKey +"_velocity";
  //const auto normalKey = boundaryKey +"_normal";
  state.enroll(pointKey,mPoint);
  state.enroll(pointKey,mVelocity);
  //state.enroll(pointKey,mNormal);
}

template<typename Dimension>
void
CylinderSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mPoint += multiplier*mVelocity;
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
}


}