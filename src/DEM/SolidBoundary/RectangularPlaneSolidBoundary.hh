//---------------------------------Spheral++----------------------------------//
// RectangularPlaneSolidBoundary -- solid planar boundary for DEM with finite
//                                  extent and rectangular shape.
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#ifndef __Spheral_RectangularPlaneSolidBoundary_hh__
#define __Spheral_RectangularPlaneSolidBoundary_hh__

#include "DEM/SolidBoundary/SolidBoundaryBase.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;

template<typename Dimension>
class RectangularPlaneSolidBoundary : public SolidBoundaryBase<Dimension> {

    typedef typename Dimension::Scalar Scalar;
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::Tensor Tensor;
    typedef typename DEMDimension<Dimension>::AngularVector RotationType;

public:
  //--------------------------- Public Interface ---------------------------//

  RectangularPlaneSolidBoundary(const Vector& point,
                                const Vector& exent, 
                                const Tensor& basis,
                                const RotationType& angularVelocity);

  ~RectangularPlaneSolidBoundary();

  virtual Vector distance(const Vector& position) const override;
  virtual Vector localVelocity(const Vector& position) const override;

  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  virtual void update(const double multiplier,
                      const double time,
                      const double dt) override;

  const Vector& point() const;
  void point(const Vector& value);

  const Tensor& basis() const;
  void basis(const Tensor& value);

  const Vector& extent() const;
  void extent(const Vector& value);

  const Vector& velocity() const;
  void velocity(const Vector& value);

  const RotationType& angularVelocity() const;
  void angularVelocity(const RotationType& value);

  virtual std::string label() const override { return "RectangularPlaneSolidBoundary" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;

protected:
  //-------------------------- Protected Interface --------------------------//
  Vector mPoint;
  Tensor mBasis;
  Vector mExtent;
  Vector mVelocity;
  RotationType mAngularVelocity;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  RectangularPlaneSolidBoundary();
  RectangularPlaneSolidBoundary(const RectangularPlaneSolidBoundary&);
  RectangularPlaneSolidBoundary& operator=(const RectangularPlaneSolidBoundary&);
};

}

#include "RectangularPlaneSolidBoundaryInline.hh"

#endif
