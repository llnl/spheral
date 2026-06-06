namespace Spheral {
//------------------------------------------------------------------------------
// The internal state field lists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const FieldList<Dimension, typename Dimension::Scalar>&
MFM<Dimension>::
DvolumeDt() const {
  return mDvolumeDt;
}
}
