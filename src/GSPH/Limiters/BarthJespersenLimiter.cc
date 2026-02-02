//---------------------------------Spheral++----------------------------------//
// BarthJespersenLimiter 
//     J. Barth, D. C. Jespersen, The design and application of upwind schemes 
//     on unstructured meshes, in: 27th Aerospace Sciences Meetings, AIAA Paper 
//     89-0366, Reno, NV, 1989. doi:10.2514/6.1989-366
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "BarthJespersenLimiter.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
BarthJespersenLimiter<Dimension>::
BarthJespersenLimiter():
  LimiterBase<Dimension>(true,true){
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
BarthJespersenLimiter<Dimension>::
~BarthJespersenLimiter(){}

//------------------------------------------------------------------------------
// slope limiter
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
BarthJespersenLimiter<Dimension>::
fluxLimiter(const typename Dimension::Scalar x) const {
  return std::min(std::min(0.5*(x+1),2.0),2.0*x) ;
}

}