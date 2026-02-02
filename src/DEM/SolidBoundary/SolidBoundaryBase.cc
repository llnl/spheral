//---------------------------------Spheral++----------------------------------//
// SolidBoundaryBase -- A base class for solid wall contact boundaries in DEM.
//                      Solid boundaries are analytically defined surfaces
//                      that behave like soft-sphere contacts. This type of
//                      boundary is distinct from standard spheral boundaries
//                      in that there are no ghost particles. The evolution of 
//                      particle-boundary contact state (i.e sliding/rolling/
//                      torsion springs) is handled exactly like 
//                      particle-particle contacts and the machinery for this
//                      is implemented in the DEMBase class.  The boundaries
//                      themselves simply define the surface and how it evolves
//                      in time.                  
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//


#include "DEM/SolidBoundary/SolidBoundaryBase.hh"

namespace Spheral {
template<typename Dimension>
SolidBoundaryBase<Dimension>::
SolidBoundaryBase():
  mUniqueIndex(-1),
  mRestart(registerWithRestart(*this)){
} 
template<typename Dimension>
SolidBoundaryBase<Dimension>::
~SolidBoundaryBase(){} 
}