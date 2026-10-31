namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
CylinderSolidBoundary<Dimension>::
point() const {
  return mPoint;
}

template<typename Dimension>
inline
void
CylinderSolidBoundary<Dimension>::
point(const typename Dimension::Vector& value) {
  mPoint=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
CylinderSolidBoundary<Dimension>::
axis() const {
  return mAxis;
}

template<typename Dimension>
inline
void
CylinderSolidBoundary<Dimension>::
axis(const typename Dimension::Vector& value) {
  mAxis=value;
}


template<typename Dimension>
inline
typename Dimension::Scalar
CylinderSolidBoundary<Dimension>::
length() const {
  return mLength;
}

template<typename Dimension>
inline
void
CylinderSolidBoundary<Dimension>::
length(typename Dimension::Scalar value) {
  mLength=value;
}


template<typename Dimension>
inline
typename Dimension::Scalar
CylinderSolidBoundary<Dimension>::
radius() const {
  return mRadius;
}

template<typename Dimension>
inline
void
CylinderSolidBoundary<Dimension>::
radius(typename Dimension::Scalar value) {
  mRadius=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
CylinderSolidBoundary<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
CylinderSolidBoundary<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}


template<typename Dimension>
inline
const typename DEMDimension<Dimension>::AngularVector&
CylinderSolidBoundary<Dimension>::
angularVelocity() const {
  return mAngularVelocity;
}

template<typename Dimension>
inline
void
CylinderSolidBoundary<Dimension>::
angularVelocity(const typename DEMDimension<Dimension>::AngularVector& value)  {
  mAngularVelocity=value;
}
}