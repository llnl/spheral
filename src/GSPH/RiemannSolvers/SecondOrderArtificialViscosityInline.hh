namespace Spheral {

//------------------------------------------------------------------------------
// set/get linear reconstruction switch
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
SecondOrderArtificialViscosity<Dimension>::
Cl() const {
  return mCl;
}

template<typename Dimension>
inline
void
SecondOrderArtificialViscosity<Dimension>::
Cl(typename Dimension::Scalar x) {
  mCl=x;
}


template<typename Dimension>
inline
typename Dimension::Scalar
SecondOrderArtificialViscosity<Dimension>::
Cq() const {
  return mCq;
}

template<typename Dimension>
inline
void
SecondOrderArtificialViscosity<Dimension>::
Cq(typename Dimension::Scalar x) {
  mCq=x;
}


}