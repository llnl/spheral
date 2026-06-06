//---------------------------------Spheral++----------------------------------//
// ClippedSphereSolidBoundary -- cylinder with finite length solid boundary for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_ClippedSphereSolidBoundary_hh__
#define __Spheral_ClippedSphereSolidBoundary_hh__

#include "DEM/SolidBoundary/SolidBoundaryBase.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class ClippedSphereSolidBoundary : public SolidBoundaryBase<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::Tensor Tensor;

public:
  //--------------------------- Public Interface ---------------------------//

  ClippedSphereSolidBoundary(const Vector& center, 
                      const Scalar  radius,
                      const Vector& clipPoint,
                      const Vector& clipAxis);

  ~ClippedSphereSolidBoundary();

  virtual Vector distance(const Vector& position) const override;
  virtual Vector localVelocity(const Vector& position) const override;

  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;


  virtual void update(const double multiplier,
                      const double time,
                      const double dt) override;

  const Vector& center() const;
  void center(const Vector& value);

  Scalar radius() const;
  void radius(Scalar value);

  const Vector& clipPoint() const;
  void clipPoint(const Vector& value);

  const Vector& clipAxis() const;
  void clipAxis(const Vector& value);

  const Vector& velocity() const;
  void velocity(const Vector& value);

  void setClipIntersectionRadius();

  virtual std::string label() const override { return "ClippedSphereSolidBoundary" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
protected:
  //-------------------------- Protected Interface --------------------------//
  Vector mCenter;
  Scalar mRadius;
  Vector mClipPoint;
  Vector mClipAxis;
  Scalar mClipIntersectionRadius;
  
  Vector mVelocity;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  ClippedSphereSolidBoundary();
  ClippedSphereSolidBoundary(const ClippedSphereSolidBoundary&);
  ClippedSphereSolidBoundary& operator=(const ClippedSphereSolidBoundary&);
};

}

#include "ClippedSphereSolidBoundaryInline.hh"

#endif
