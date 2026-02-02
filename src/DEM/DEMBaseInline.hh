namespace Spheral {

//------------------------------------------------------------------------------
// Access the optional min & max bounds for generating meshes.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const typename Dimension::Vector&
DEMBase<Dimension>::
xmin() const {
  return mxmin;
}

template<typename Dimension>
inline
const typename Dimension::Vector&
DEMBase<Dimension>::
xmax() const {
  return mxmax;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
xmin(const typename Dimension::Vector& x) {
  mxmin = x;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
xmax(const typename Dimension::Vector& x) {
  mxmax = x;
}


//------------------------------------------------------------------------------
// set get for numerical switches
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
DEMBase<Dimension>::
contactRemovalFrequency() const {
  return mContactRemovalFrequency;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
contactRemovalFrequency(int x) {
  mContactRemovalFrequency = x;
}

template<typename Dimension>
inline
int
DEMBase<Dimension>::
cycle() const {
  return mCycle;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
cycle(int x) {
  mCycle = x;
}


//------------------------------------------------------------------------------
// CFL number (ratio to estimated contact duration)
//------------------------------------------------------------------------------

template<typename Dimension>
inline
typename Dimension::Scalar
DEMBase<Dimension>::
stepsPerCollision() const {
  return mStepsPerCollision;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
stepsPerCollision(typename Dimension::Scalar x) {
  mStepsPerCollision = x;
}

//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, int>&
DEMBase<Dimension>::
timeStepMask() const {
  return mTimeStepMask;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
DEMBase<Dimension>::
DxDt() const {
  return mDxDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Vector>&
DEMBase<Dimension>::
DvDt() const {
  return mDvDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename DEMDimension<Dimension>::AngularVector>&
DEMBase<Dimension>::
DomegaDt() const {
  return mDomegaDt;
}

template<typename Dimension>
inline
const FieldList<Dimension, typename DEMDimension<Dimension>::AngularVector>&
DEMBase<Dimension>::
omega() const {
  return mOmega;
}


//------------------------------------------------------------------------------
// Pair things
//------------------------------------------------------------------------------


template<typename Dimension>
inline
const FieldList<Dimension, vector<int>>&
DEMBase<Dimension>::
isActiveContact() const {
  return mIsActiveContact;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<size_t>>&
DEMBase<Dimension>::
neighborIndices() const {
  return mNeighborIndices;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
shearDisplacement() const {
  return mShearDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
rollingDisplacement() const {
  return mRollingDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Scalar>>&
DEMBase<Dimension>::
torsionalDisplacement() const {
  return mTorsionalDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
DDtShearDisplacement() const {
  return mDDtShearDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
newShearDisplacement() const {
  return mNewShearDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
DDtRollingDisplacement() const {
  return mDDtRollingDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Vector>>&
DEMBase<Dimension>::
newRollingDisplacement() const {
  return mNewRollingDisplacement;
}


template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Scalar>>&
DEMBase<Dimension>::
DDtTorsionalDisplacement() const {
  return mDDtTorsionalDisplacement;
}

template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Scalar>>&
DEMBase<Dimension>::
newTorsionalDisplacement() const {
  return mNewTorsionalDisplacement;
}


template<typename Dimension>
inline
const FieldList<Dimension, vector<typename Dimension::Scalar>>&
DEMBase<Dimension>::
equilibriumOverlap() const {
  return mEquilibriumOverlap;
}

template<typename Dimension>
inline
const vector<ContactIndex>&
DEMBase<Dimension>::
contactStorageIndices() const {
  return mContactStorageIndices;
}

template<typename Dimension>
inline
const DataBase<Dimension>&
DEMBase<Dimension>::
dataBase() const {
  return mDataBase;
}
//------------------------------------------------------------------------------
// torsion specializations
//------------------------------------------------------------------------------
template<>
inline
DEMDimension<Dim<1>>::AngularVector
DEMBase<Dim<1>>::
torsionMoment(const Dim<1>::Vector rhatij, 
              const DEMDimension<Dim<1>>::AngularVector omegai,
              const DEMDimension<Dim<1>>::AngularVector omegaj) const {
  return 0.0;
}

template<>
inline
DEMDimension<Dim<2>>::AngularVector
DEMBase<Dim<2>>::
torsionMoment(const Dim<2>::Vector rhatij, 
              const DEMDimension<Dim<2>>::AngularVector omegai,
              const DEMDimension<Dim<2>>::AngularVector omegaj) const {
  return 0.0;
}

template<>
inline
DEMDimension<Dim<3>>::AngularVector
DEMBase<Dim<3>>::
torsionMoment(const Dim<3>::Vector rhatij, 
              const DEMDimension<Dim<3>>::AngularVector omegai,
              const DEMDimension<Dim<3>>::AngularVector omegaj) const {
  return rhatij;
}


//------------------------------------------------------------------------------
// rolling Moment specializations
//------------------------------------------------------------------------------
template<>
inline
DEMDimension<Dim<1>>::AngularVector
DEMBase<Dim<1>>::
rollingMoment(const Dim<1>::Vector rhatij, 
              const Dim<1>::Vector vroti,
              const Dim<1>::Vector vrotj) const {
  return 0.0;
}

template<>
inline
DEMDimension<Dim<2>>::AngularVector
DEMBase<Dim<2>>::
rollingMoment(const Dim<2>::Vector rhatij, 
              const Dim<2>::Vector vroti,
              const Dim<2>::Vector vrotj) const {
  return ( DEMDimension<Dim<2>>::cross((vroti + vrotj),rhatij) > 0.0 ? 1.0 : -1.0 );
}

template<>
inline
DEMDimension<Dim<3>>::AngularVector
DEMBase<Dim<3>>::
rollingMoment(const Dim<3>::Vector rhatij, 
              const Dim<3>::Vector vroti,
              const Dim<3>::Vector vrotj) const {
  return DEMDimension<Dim<3>>::cross((vroti + vrotj),rhatij).unitVector();
}

//------------------------------------------------------------------------------
// Add a Boundary condition to the end of the current boundary list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
DEMBase<Dimension>::
appendSolidBoundary(SolidBoundaryBase<Dimension>& boundary) {
  mNewSolidBoundaryIndex -= 1;
  boundary.uniqueIndex(mNewSolidBoundaryIndex);
  mSolidBoundaries.push_back(&boundary);
}

//------------------------------------------------------------------------------
// Clear (erase) the boundary condition list.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
DEMBase<Dimension>::
clearSolidBoundaries() {
  mSolidBoundaries = std::vector<SolidBoundaryBase<Dimension>*>();
}


template<typename Dimension>
inline
int
DEMBase<Dimension>::
newSolidBoundaryIndex() const {
  return mNewSolidBoundaryIndex;
}

//------------------------------------------------------------------------------
// Test if the given Boundary condition is listed in the physics package.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
DEMBase<Dimension>::
haveSolidBoundary(const SolidBoundaryBase<Dimension>& boundary) const {
  return std::count(mSolidBoundaries.begin(), mSolidBoundaries.end(), &boundary) > 0;
}

template<typename Dimension>
inline
void
DEMBase<Dimension>::
removeSolidBoundary(const SolidBoundaryBase<Dimension>& boundary) {
  const auto bcPtr = std::find(mSolidBoundaries.begin(),mSolidBoundaries.end(),&boundary);
  if (bcPtr != mSolidBoundaries.end()) mSolidBoundaries.erase(bcPtr);
}


template<typename Dimension>
inline
unsigned int
DEMBase<Dimension>::
numSolidBoundaries() const {
  return mSolidBoundaries.size();
}

template<typename Dimension>
inline
const std::vector<SolidBoundaryBase<Dimension>*>&
DEMBase<Dimension>::solidBoundaryConditions() const {
  return mSolidBoundaries;
}


//------------------------------------------------------------------------------
// contact counts
//------------------------------------------------------------------------------
template<typename Dimension>
inline
unsigned int
DEMBase<Dimension>::
numParticleParticleContacts() const {
  const auto& cm = mDataBase.connectivityMap();
  const auto& pairs = cm.nodePairList();
  return pairs.size();
}

template<typename Dimension>
inline
unsigned int
DEMBase<Dimension>::
numContacts() const {
  return mContactStorageIndices.size();
}

template<typename Dimension>
inline
unsigned int
DEMBase<Dimension>::
numParticleBoundaryContacts() const {
  return (mContactStorageIndices.size()-this->numParticleParticleContacts());
}

}
