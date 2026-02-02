namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
SphereSolidBoundary<Dimension>::
center() const {
  return mCenter;
}

template<typename Dimension>
inline
void
SphereSolidBoundary<Dimension>::
center(const typename Dimension::Vector& value) {
  mCenter=value;
}

template<typename Dimension>
inline
typename Dimension::Scalar
SphereSolidBoundary<Dimension>::
radius() const {
  return mRadius;
}

template<typename Dimension>
inline
void
SphereSolidBoundary<Dimension>::
radius(typename Dimension::Scalar value) {
  mRadius=value;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
SphereSolidBoundary<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
SphereSolidBoundary<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}

template<typename Dimension>
inline
const typename DEMDimension<Dimension>::AngularVector&
SphereSolidBoundary<Dimension>::
angularVelocity() const {
  return mAngularVelocity;
}

template<typename Dimension>
inline
void
SphereSolidBoundary<Dimension>::
angularVelocity(const typename DEMDimension<Dimension>::AngularVector& value)  {
  mAngularVelocity=value;
}

}