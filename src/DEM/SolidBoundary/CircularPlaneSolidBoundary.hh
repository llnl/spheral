//---------------------------------Spheral++----------------------------------//
// CircularPlaneSolidBoundary -- solid planar boundary for DEM with finite 
//                               extent and circular shape.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_CircularPlaneSolidBoundary_hh__
#define __Spheral_CircularPlaneSolidBoundary_hh__

#include "DEM/SolidBoundary/SolidBoundaryBase.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class CircularPlaneSolidBoundary : public SolidBoundaryBase<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;
    typedef typename DEMDimension<Dimension>::AngularVector RotationType;

public:
  //--------------------------- Public Interface ---------------------------//

  CircularPlaneSolidBoundary(const Vector& point,
                             const Vector& normal,
                             const Scalar& exent,
                             const RotationType& angularVelocity);

  ~CircularPlaneSolidBoundary();

  virtual Vector distance(const Vector& position) const override;
  virtual Vector localVelocity(const Vector& position) const override;

  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  virtual void update(const double multiplier,
                      const double time,
                      const double dt) override;

  const Vector& point() const;
  void point(const Vector& value);

    const Vector& normal() const;
  void normal(const Vector& value);

  Scalar extent() const;
  void extent(Scalar value);

  const Vector& velocity() const;
  void velocity(const Vector& value);

  const RotationType& angularVelocity() const;
  void angularVelocity(const RotationType& value);

  virtual std::string label() const override { return "CircularPlaneSolidBoundary" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;

protected:
  //-------------------------- Protected Interface --------------------------//
  Vector mPoint;
  Vector mNormal;
  Scalar mExtent;
  Vector mVelocity;
  RotationType mAngularVelocity;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  CircularPlaneSolidBoundary();
  CircularPlaneSolidBoundary(const CircularPlaneSolidBoundary&);
  CircularPlaneSolidBoundary& operator=(const CircularPlaneSolidBoundary&);
};

}

#include "CircularPlaneSolidBoundaryInline.hh"

#endif
