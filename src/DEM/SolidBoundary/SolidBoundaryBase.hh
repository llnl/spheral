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
//----------------------------------------------------------------------------//
// ToDo
// -- add complete registration methods
// -- set the unique index of the bc so no std::str input req. for regState
// -- make restartable 
// -- stateBase add scalar and tensor (talk to Mike about different pattern)
// -- tests
//-----------------------------------------------------------------------------
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_SolidBoundaryBase_hh__
#define __Spheral_SolidBoundaryBase_hh__

#include "DataOutput/registerWithRestart.hh"

#include <string>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
class FileIO;

template<typename Dimension>
class SolidBoundaryBase {

public:
//--------------------------- Public Interface ---------------------------//

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  SolidBoundaryBase();

  virtual ~SolidBoundaryBase();

  virtual Vector distance(const Vector& position) const = 0;
  virtual Vector localVelocity(const Vector& position) const = 0;

  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) = 0;

  virtual void update(const double multiplier,
                      const double t,
                      const double dt) = 0;
  
  void uniqueIndex(int uId);
  int uniqueIndex() const;

  // restartability will default to no-op 
  virtual std::string label() const { return "SolidBoundaryBase" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const {};
  virtual void restoreState(const FileIO& file, const std::string& pathName) {};

private:
//--------------------------- Public Interface ---------------------------//
int mUniqueIndex;

RestartRegistrationType mRestart;

};
}

#include "SolidBoundaryBaseInline.hh"

#endif
