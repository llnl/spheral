namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
ClippedSphereSolidBoundary<Dimension>::
center() const {
  return mCenter;
}

template<typename Dimension>
inline
void
ClippedSphereSolidBoundary<Dimension>::
center(const typename Dimension::Vector& value) {
  mCenter=value;
}

template<typename Dimension>
inline
typename Dimension::Scalar
ClippedSphereSolidBoundary<Dimension>::
radius() const {
  return mRadius;
}

template<typename Dimension>
inline
void
ClippedSphereSolidBoundary<Dimension>::
radius(typename Dimension::Scalar value) {
  mRadius=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
ClippedSphereSolidBoundary<Dimension>::
clipPoint() const {
  return mClipPoint;
}

template<typename Dimension>
inline
void
ClippedSphereSolidBoundary<Dimension>::
clipPoint(const typename Dimension::Vector& value) {
  mClipPoint=value;
  this->setClipIntersectionRadius();
}


template<typename Dimension>
inline
const typename Dimension::Vector&
ClippedSphereSolidBoundary<Dimension>::
clipAxis() const {
  return mClipAxis;
}

template<typename Dimension>
inline
void
ClippedSphereSolidBoundary<Dimension>::
clipAxis(const typename Dimension::Vector& value) {
  mClipAxis=value.unitVector();
  this->setClipIntersectionRadius();
}

template<typename Dimension>
inline
const typename Dimension::Vector&
ClippedSphereSolidBoundary<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
ClippedSphereSolidBoundary<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}

}