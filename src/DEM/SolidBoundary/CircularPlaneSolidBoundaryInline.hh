namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
CircularPlaneSolidBoundary<Dimension>::
point() const {
  return mPoint;
}

template<typename Dimension>
inline
void
CircularPlaneSolidBoundary<Dimension>::
point(const typename Dimension::Vector& value) {
  mPoint=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
CircularPlaneSolidBoundary<Dimension>::
normal() const {
  return mNormal;
}

template<typename Dimension>
inline
void
CircularPlaneSolidBoundary<Dimension>::
normal(const typename Dimension::Vector& value) {
  mNormal = value;
}


template<typename Dimension>
inline
typename Dimension::Scalar
CircularPlaneSolidBoundary<Dimension>::
extent() const {
  return mExtent;
}

template<typename Dimension>
inline
void
CircularPlaneSolidBoundary<Dimension>::
extent(typename Dimension::Scalar value) {
  mExtent=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
CircularPlaneSolidBoundary<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
CircularPlaneSolidBoundary<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}


template<typename Dimension>
inline
const typename DEMDimension<Dimension>::AngularVector&
CircularPlaneSolidBoundary<Dimension>::
angularVelocity() const {
  return mAngularVelocity;
}

template<typename Dimension>
inline
void
CircularPlaneSolidBoundary<Dimension>::
angularVelocity(const typename DEMDimension<Dimension>::AngularVector& value) {
  mAngularVelocity=value;
}

}