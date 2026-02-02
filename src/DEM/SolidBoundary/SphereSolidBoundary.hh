//---------------------------------Spheral++----------------------------------//
// SphereSolidBoundary -- cylinder with finite length solid boundary for DEM
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_SphereSolidBoundary_hh__
#define __Spheral_SphereSolidBoundary_hh__

#include "DEM/DEMDimension.hh"
#include "DEM/SolidBoundary/SolidBoundaryBase.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class SphereSolidBoundary : public SolidBoundaryBase<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::Tensor Tensor;
    typedef typename DEMDimension<Dimension>::AngularVector RotationType;
public:
  //--------------------------- Public Interface ---------------------------//
  SphereSolidBoundary(const Vector& center, 
                      const Scalar  radius,
                      const RotationType& angularVelocity);

  ~SphereSolidBoundary();

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

  const Vector& velocity() const;
  void velocity(const Vector& value);

  const RotationType& angularVelocity() const;
  void angularVelocity(const RotationType& value);

  virtual std::string label() const override { return "SphereSolidBoundary" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;

protected:
  //-------------------------- Protected Interface --------------------------//
  Vector mCenter;
  Scalar mRadius;
  
  Vector mVelocity;

  RotationType mAngularVelocity;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  SphereSolidBoundary();
  SphereSolidBoundary(const SphereSolidBoundary&);
  SphereSolidBoundary& operator=(const SphereSolidBoundary&);
};

}

#include "SphereSolidBoundaryInline.hh"

#endif
