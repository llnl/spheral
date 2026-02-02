//---------------------------------Spheral++----------------------------------//
// BarthJespersenLimiter 
//     J. Barth, D. C. Jespersen, The design and application of upwind schemes 
//     on unstructured meshes, in: 27th Aerospace Sciences Meetings, AIAA Paper 
//     89-0366, Reno, NV, 1989. doi:10.2514/6.1989-366
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_BarthJespersenLimiter_hh__
#define __Spheral_BarthJespersenLimiter_hh__

#include "LimiterBase.hh"

namespace Spheral {

template<typename Dimension>
class BarthJespersenLimiter : public LimiterBase<Dimension> {

public:

  typedef typename Dimension::Scalar Scalar;

  BarthJespersenLimiter();

  ~BarthJespersenLimiter();

  virtual
  Scalar fluxLimiter(const Scalar) const override;

};


}

#else

// Forward declaration.
namespace Spheral {
  template<typename Dimension> class BarthJespersenLimiter;
}

#endif

