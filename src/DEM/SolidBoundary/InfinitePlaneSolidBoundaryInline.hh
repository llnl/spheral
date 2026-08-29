namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
InfinitePlaneSolidBoundary<Dimension>::
point() const {
  return mPoint;
}

template<typename Dimension>
inline
void
InfinitePlaneSolidBoundary<Dimension>::
point(const typename Dimension::Vector& value) {
  mPoint=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
InfinitePlaneSolidBoundary<Dimension>::
normal() const {
  return mNormal;
}

template<typename Dimension>
inline
void
InfinitePlaneSolidBoundary<Dimension>::
normal(const typename Dimension::Vector& value) {
  mNormal = value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
InfinitePlaneSolidBoundary<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
InfinitePlaneSolidBoundary<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}


template<typename Dimension>
inline
const typename DEMDimension<Dimension>::AngularVector&
InfinitePlaneSolidBoundary<Dimension>::
angularVelocity() const {
  return mAngularVelocity;
}

template<typename Dimension>
inline
void
InfinitePlaneSolidBoundary<Dimension>::
angularVelocity(const typename DEMDimension<Dimension>::AngularVector& value)  {
  mAngularVelocity=value;
}

}