namespace Spheral {

template<typename Dimension>
inline
const typename Dimension::Vector&
RectangularPlaneSolidBoundary<Dimension>::
point() const {
  return mPoint;
}

template<typename Dimension>
inline
void
RectangularPlaneSolidBoundary<Dimension>::
point(const typename Dimension::Vector& value) {
  mPoint=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
RectangularPlaneSolidBoundary<Dimension>::
extent() const {
  return mExtent;
}

template<typename Dimension>
inline
void
RectangularPlaneSolidBoundary<Dimension>::
extent(const typename Dimension::Vector& value) {
  mExtent=value;
}


template<typename Dimension>
inline
const typename Dimension::Tensor&
RectangularPlaneSolidBoundary<Dimension>::
basis() const {
  return mBasis;
}

template<typename Dimension>
inline
void
RectangularPlaneSolidBoundary<Dimension>::
basis(const typename Dimension::Tensor& value) {
  mBasis=value;
}


template<typename Dimension>
inline
const typename Dimension::Vector&
RectangularPlaneSolidBoundary<Dimension>::
velocity() const {
  return mVelocity;
}

template<typename Dimension>
inline
void
RectangularPlaneSolidBoundary<Dimension>::
velocity(const typename Dimension::Vector& value)  {
  mVelocity=value;
}


template<typename Dimension>
inline
const typename DEMDimension<Dimension>::AngularVector&
RectangularPlaneSolidBoundary<Dimension>::
angularVelocity() const {
  return mAngularVelocity;
}

template<typename Dimension>
inline
void
RectangularPlaneSolidBoundary<Dimension>::
angularVelocity(const typename DEMDimension<Dimension>::AngularVector& value)  {
  mAngularVelocity=value;
}

}