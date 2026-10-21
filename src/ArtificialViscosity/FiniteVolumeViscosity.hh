//---------------------------------Spheral++----------------------------------//
// A finite-volume based viscosity.  Assumes you have constructed the
// tessellation in the state.
//
// Created by JMO, Tue Aug 13 09:43:37 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_FiniteVolumeViscosity__
#define __Spheral_FiniteVolumeViscosity__

#include "ArtificialViscosity.hh"
#include "ArtificialViscosityVariant.hh"

namespace Spheral {

template<typename Dimension>
class FiniteVolumeViscosity: public ArtificialViscosity<Dimension>,
                             public FiniteVolumeViscosityView<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ViewType = FiniteVolumeViscosityView<Dimension>;

  // Constructor, destructor
  FiniteVolumeViscosity(const Scalar Clinear,
                        const Scalar Cquadratic,
                        const TableKernel<Dimension>& WT) :
    ArtificialViscosity<Dimension>(WT),
    ViewType(Clinear, Cquadratic) { }

  virtual ~FiniteVolumeViscosity() = default;

  // No default construction, copying, or assignment
  FiniteVolumeViscosity() = delete;
  FiniteVolumeViscosity(const FiniteVolumeViscosity&) = delete;
  FiniteVolumeViscosity& operator=(const FiniteVolumeViscosity&) const = delete;

  // We are going to use a velocity gradient
  virtual bool requireVelocityGradient() const override { return true; }

  // Restart methods.
  virtual std::string label()            const override { return "FiniteVolumeViscosity"; }

  // Override the method of computing the velocity gradient
  virtual void updateVelocityGradient(const DataBase<Dimension>& db,
                                      const State<Dimension>& state,
                                      const StateDerivatives<Dimension>& derivs) override;

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
