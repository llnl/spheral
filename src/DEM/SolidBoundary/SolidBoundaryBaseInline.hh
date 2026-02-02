namespace Spheral {
//------------------------------------------------------------------------------
// set/get unique index for the solid bc
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
SolidBoundaryBase<Dimension>::
uniqueIndex(int uId){
  mUniqueIndex = uId;
}
template<typename Dimension>
inline
int
SolidBoundaryBase<Dimension>::
uniqueIndex() const {
  return mUniqueIndex;
}
}

