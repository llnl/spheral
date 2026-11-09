//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_LimitedMonaghanGingoldViscosity__
#define __Spheral_LimitedMonaghanGingoldViscosity__

#include "ArtificialViscosity.hh"
#include "ArtificialViscosityVariant.hh"

namespace Spheral {

template<typename Dimension>
class LimitedMonaghanGingoldViscosity: public ArtificialViscosity<Dimension>,
                                       public LimitedMonaghanGingoldViscosityView<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ViewType = LimitedMonaghanGingoldViscosityView<Dimension>;

  // Constructors.
  LimitedMonaghanGingoldViscosity(const Scalar Clinear,
                                  const Scalar Cquadratic,
                                  const TableKernel<Dimension>& kernel,
                                  const bool linearInExpansion,
                                  const bool quadraticInExpansion,
                                  const Scalar etaCritFrac,
                                  const Scalar etaFoldFrac) :
    ArtificialViscosity<Dimension>(kernel),
    ViewType(Clinear, Cquadratic,
             linearInExpansion, quadraticInExpansion,
             etaCritFrac, etaFoldFrac) { }

  virtual ~LimitedMonaghanGingoldViscosity() = default;

  // No default construction, copying, or assignment
  LimitedMonaghanGingoldViscosity() = delete;
  LimitedMonaghanGingoldViscosity(const LimitedMonaghanGingoldViscosity&) = delete;
  LimitedMonaghanGingoldViscosity& operator=(const LimitedMonaghanGingoldViscosity&) = delete;

  // We need the velocity gradient
  virtual bool requireVelocityGradient() const override { return true; }

  // Restart methods.
  virtual std::string label() const override { return "LimitedMonaghanGingoldViscosity"; }

  // Forward ArtificialViscosity's virtual parameter interface to ViewType.
  virtual Scalar Cl() const override { return ViewType::Cl(); }
  virtual Scalar Cq() const override { return ViewType::Cq(); }
  virtual bool balsaraShearCorrection() const override {
    return ViewType::balsaraShearCorrection();
  }
  virtual Scalar epsilon2() const override { return ViewType::epsilon2(); }
  virtual Scalar negligibleSoundSpeed() const override {
    return ViewType::negligibleSoundSpeed();
  }

  virtual void Cl(const Scalar x) override { ViewType::Cl(x); }
  virtual void Cq(const Scalar x) override { ViewType::Cq(x); }
  virtual void balsaraShearCorrection(const bool x) override {
    ViewType::balsaraShearCorrection(x);
  }
  virtual void epsilon2(const Scalar x) override { ViewType::epsilon2(x); }
  virtual void negligibleSoundSpeed(const Scalar x) override {
    ViewType::negligibleSoundSpeed(x);
  }

  Scalar etaCritFrac() const { return ViewType::etaCritFrac(); }
  Scalar etaFoldFrac() const { return ViewType::etaFoldFrac(); }
  void etaCritFrac(const Scalar x) { ViewType::etaCritFrac(x); }
  void etaFoldFrac(const Scalar x) { ViewType::etaFoldFrac(x); }

  // Return a device-safe snapshot of the inherited view state.
  ViewType view() const {
    return static_cast<const ViewType&>(*this);
  }
  ArtificialViscosityVariant<Dimension> variantView() const override {
    return this->view();
  }
};

}

#endif
