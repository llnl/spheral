//---------------------------------Spheral++----------------------------------//
// ClippedSphereSolidBoundary -- N-dimensional spherical solid boundary for DEM.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "FileIO/FileIO.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"

#include "DEM/SolidBoundary/ClippedSphereSolidBoundary.hh"

#include <cmath>

#include <string>
using std::string;

namespace Spheral {

template<typename Dimension>
ClippedSphereSolidBoundary<Dimension>::
ClippedSphereSolidBoundary(const Vector& center,
                    const Scalar  radius,
                    const Vector& clipPoint,
                    const Vector& clipAxis):
  SolidBoundaryBase<Dimension>(),
  mCenter(center),
  mRadius(radius),
  mClipPoint(clipPoint),
  mClipAxis(clipAxis),
  mClipIntersectionRadius(0.0),
  mVelocity(Vector::zero()){
    this->setClipIntersectionRadius();
    mClipAxis = mClipAxis.unitVector();
}

template<typename Dimension>
ClippedSphereSolidBoundary<Dimension>::
~ClippedSphereSolidBoundary(){
}


template<typename Dimension>
typename Dimension::Vector
ClippedSphereSolidBoundary<Dimension>::
distance(const Vector& position) const { 

  // contacting sphere
  const auto contactPoint = (position - mCenter).unitVector()*mRadius + mCenter;
  Vector dist = position - contactPoint;

  const auto planeSignedDistance = (contactPoint - mClipPoint).dot(mClipAxis);

  // if contant pt above clip plane check for edge contact
  if (planeSignedDistance > 0.0){

    // break into perp and in-plane components
    const auto q = (position-mClipPoint);
    const auto qnMag =  q.dot(mClipAxis);
    const auto qn = qnMag * mClipAxis;
    const auto qr = q - qn;

    // if outside circle enforce planar solid bc
    dist = min(qr.magnitude() - mClipIntersectionRadius,0.0)*qr.unitVector() + qn;

  }
  return dist;
}

template<typename Dimension>
typename Dimension::Vector
ClippedSphereSolidBoundary<Dimension>::
localVelocity(const Vector& position) const { 
  return mVelocity;
}

template<typename Dimension>
void
ClippedSphereSolidBoundary<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  const auto boundaryKey = "ClippedSphereSolidBoundary_" + std::to_string(std::abs(this->uniqueIndex()));
  const auto pointKey = boundaryKey +"_point";
  const auto clipPointKey = boundaryKey +"_clipPoint";
  const auto velocityKey = boundaryKey +"_velocity";

  state.enroll(pointKey,mCenter);
  state.enroll(clipPointKey,mClipPoint);
  state.enroll(pointKey,mVelocity);

}

template<typename Dimension>
void
ClippedSphereSolidBoundary<Dimension>::
update(const double multiplier, const double t, const double dt) {   
  mCenter += multiplier*mVelocity;
  mClipPoint += multiplier*mVelocity;
}


template<typename Dimension>
void
ClippedSphereSolidBoundary<Dimension>::
setClipIntersectionRadius() {
  const auto rcMag = (mClipPoint - mCenter).dot(mClipAxis);
  mClipIntersectionRadius = (rcMag < mRadius ? std::sqrt(mRadius*mRadius-rcMag*rcMag) : 0.0);
  mClipPoint = rcMag * mClipAxis + mCenter;
}

//------------------------------------------------------------------------------
// Restart
//------------------------------------------------------------------------------
template<typename Dimension>
void
ClippedSphereSolidBoundary<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mCenter, pathName + "/center");
  file.write(mRadius, pathName + "/radius");
  file.write(mClipPoint, pathName + "/clipPoint");
  file.write(mClipAxis, pathName + "/clipAxis");
  file.write(mClipIntersectionRadius, pathName + "/clipIntersectionRadius");
  file.write(mVelocity, pathName + "/velocity");
}


template<typename Dimension>
void
ClippedSphereSolidBoundary<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mCenter, pathName + "/center");
  file.read(mRadius, pathName + "/radius");
  file.read(mClipPoint, pathName + "/clipPoint");
  file.read(mClipAxis, pathName + "/clipAxis");
  file.read(mClipIntersectionRadius, pathName + "/clipIntersectionRadius");
  file.read(mVelocity, pathName + "/velocity");
}

}