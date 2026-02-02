namespace Spheral {


//------------------------------------------------------------------------------
// set/get to activate/deactivate fast timestepping
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
LinearSpringDEM<Dimension>::
enableFastTimeStepping() const {
  return mEnableFastTimeStepping;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
enableFastTimeStepping(bool x) {
  mEnableFastTimeStepping = x;
}


//------------------------------------------------------------------------------
// set/get our spring constant
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
normalSpringConstant() const {
  return mNormalSpringConstant;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
normalSpringConstant(typename Dimension::Scalar x) {
  mNormalSpringConstant = x;
}


//------------------------------------------------------------------------------
// set/get the restitution coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
normalRestitutionCoefficient() const {
  return mNormalRestitutionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
normalRestitutionCoefficient(typename Dimension::Scalar x) {
  mNormalRestitutionCoefficient = x;
}

//------------------------------------------------------------------------------
// set/get our spring constant
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
tangentialSpringConstant() const {
  return mTangentialSpringConstant;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
tangentialSpringConstant(typename Dimension::Scalar x) {
  mTangentialSpringConstant = x;
}


//------------------------------------------------------------------------------
// set/get the restitution coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
tangentialRestitutionCoefficient() const {
  return mTangentialRestitutionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
tangentialRestitutionCoefficient(typename Dimension::Scalar x) {
  mTangentialRestitutionCoefficient = x;
}

//------------------------------------------------------------------------------
// set/get the dynamic friction coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
dynamicFrictionCoefficient() const {
  return mDynamicFrictionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
dynamicFrictionCoefficient(typename Dimension::Scalar x) {
  mDynamicFrictionCoefficient = x;
}


//------------------------------------------------------------------------------
// set/get the static friction coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
staticFrictionCoefficient() const {
  return mStaticFrictionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
staticFrictionCoefficient(typename Dimension::Scalar x) {
  mStaticFrictionCoefficient = x;
}

//------------------------------------------------------------------------------
// set/get the rolling friction coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
rollingFrictionCoefficient() const {
  return mRollingFrictionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
rollingFrictionCoefficient(typename Dimension::Scalar x) {
  mRollingFrictionCoefficient = x;
}


//------------------------------------------------------------------------------
// set/get the torsional friction coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
torsionalFrictionCoefficient() const {
  return mTorsionalFrictionCoefficient;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
torsionalFrictionCoefficient(typename Dimension::Scalar x) {
  mTorsionalFrictionCoefficient = x;
}

//------------------------------------------------------------------------------
// cohesive coefficient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
cohesiveTensileStrength() const {
  return mCohesiveTensileStrength;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
cohesiveTensileStrength(typename Dimension::Scalar x) {
  mCohesiveTensileStrength = x;
}

//------------------------------------------------------------------------------
// set/get the shape factor
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
shapeFactor() const {
  return mShapeFactor;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
shapeFactor(typename Dimension::Scalar x) {
  mShapeFactor = x;
}

//------------------------------------------------------------------------------
// set/get this convienent damping parameter
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
normalBeta() const {
  return mNormalBeta;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
normalBeta(typename Dimension::Scalar x) {
  mNormalBeta = x;
}

template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
tangentialBeta() const {
  return mTangentialBeta;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
tangentialBeta(typename Dimension::Scalar x) {
  mTangentialBeta = x;
}

template<typename Dimension>
inline
typename Dimension::Scalar
LinearSpringDEM<Dimension>::
collisionDuration() const {
  return mCollisionDuration;
}
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
collisionDuration(typename Dimension::Scalar x) {
  mCollisionDuration = x;
}

//------------------------------------------------------------------------------
// moment of interia specializations
//------------------------------------------------------------------------------
template<>
inline
Dim<1>::Scalar
LinearSpringDEM<Dim<1>>::
momentOfInertia(const Dim<1>::Scalar m, const Dim<1>::Scalar R) const {
  return 0.5*m*R*R;
}

template<>
inline
Dim<2>::Scalar
LinearSpringDEM<Dim<2>>::
momentOfInertia(const Dim<2>::Scalar m, const Dim<2>::Scalar R) const {
  return 0.5*m*R*R;
}

template<>
inline
Dim<3>::Scalar
LinearSpringDEM<Dim<3>>::
momentOfInertia(const Dim<3>::Scalar m, const Dim<3>::Scalar R) const {
  return 0.4*m*R*R;
}


//------------------------------------------------------------------------------
// friction functions
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
slidingSpringDamper(const typename Dimension::Scalar k,
                    const typename Dimension::Scalar C,
                    const typename Dimension::Scalar mus,
                    const typename Dimension::Scalar mud,
                    const typename Dimension::Vector& x,
                    const typename Dimension::Vector& DxDt,
                    const typename Dimension::Scalar fnMag,
                    const typename Dimension::Scalar invK,
                    const typename Dimension::Vector& rhatij,
                    const bool allowSliding,
                          typename Dimension::Vector& xNew,
                          typename Dimension::Vector& forceTotal) const{

  xNew = (x - rhatij.dot(x)*rhatij).unitVector()*x.magnitude();
  
  const Vector forceSpring = - k * xNew;
  const Vector forceDamper = - C * DxDt;
  
  forceTotal  = forceSpring + forceDamper;

  if (allowSliding and (forceTotal.magnitude() > mus * fnMag)){
    forceTotal = mud*fnMag*forceTotal.unitVector();
    xNew = (forceDamper.magnitude() > mud*fnMag ? 
                      Vector::zero() : 
                      -(forceTotal-forceDamper)*invK );
  }
}

template<typename Dimension>
inline
void
LinearSpringDEM<Dimension>::
slidingSpringDamper(const typename Dimension::Scalar k,
                    const typename Dimension::Scalar C,
                    const typename Dimension::Scalar mus,
                    const typename Dimension::Scalar mud,
                    const typename Dimension::Scalar x,
                    const typename Dimension::Scalar DxDt,
                    const typename Dimension::Scalar fnMag,
                    const typename Dimension::Scalar invK,
                    const bool allowSliding,
                          typename Dimension::Scalar& xNew,
                          typename Dimension::Scalar& forceTotal) const{
  xNew = x;

  const Scalar forceSpring = - k * xNew;
  const Scalar forceDamper = - C * DxDt;
  
  forceTotal  = forceSpring + forceDamper;

  if (allowSliding and (std::abs(forceTotal) > mus * fnMag)){
    forceTotal = (forceTotal > 0.0 ? 1.0 : -1.0) * mud * fnMag;
    xNew = (std::abs(forceDamper) > mud * fnMag ? 0.0 : -(forceTotal-forceDamper)*invK);
  } 
}

//------------------------------------------------------------------------------
// FieldList
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
LinearSpringDEM<Dimension>::
momentOfInertia() const {
  return mMomentOfInertia;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
LinearSpringDEM<Dimension>::
maximumOverlap() const {
  return mMaximumOverlap;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
LinearSpringDEM<Dimension>::
newMaximumOverlap() const {
  return mNewMaximumOverlap;
}

}
