//---------------------------------Spheral++----------------------------------//
// A modified form of the Monaghan & Gingold viscosity, extended to tensor
// formalism.
//
// Created by J. Michael Owen, Mon Sep  2 14:45:35 PDT 2002
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorMonaghanGingoldViscosity__
#define __Spheral_TensorMonaghanGingoldViscosity__

#include "ArtificialViscosity.hh"
#include "TensorMonaghanGingoldViscosityView.hh"

namespace Spheral {

template<typename Dimension>
class TensorMonaghanGingoldViscosity : public ArtificialViscosity<Dimension>,
                                       public TensorMonaghanGingoldViscosityView<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ViewType = TensorMonaghanGingoldViscosityView<Dimension>;

  // Constructors and destuctor
  TensorMonaghanGingoldViscosity(const Scalar Clinear,
                                 const Scalar Cquadratic,
                                 const TableKernel<Dimension>& kernel) :
    ArtificialViscosity<Dimension>(kernel),
    ViewType(Clinear, Cquadratic) { }

  virtual ~TensorMonaghanGingoldViscosity() = default;

  // No default construction, copying, or assignment
  TensorMonaghanGingoldViscosity() = delete;
  TensorMonaghanGingoldViscosity(const TensorMonaghanGingoldViscosity&) = delete;
  TensorMonaghanGingoldViscosity& operator=(const TensorMonaghanGingoldViscosity&) = delete;

  // We need the velocity gradient
  virtual bool requireVelocityGradient() const override { return true; }

  // Restart methods.
  virtual std::string label() const override { return "TensorMonaghanGingoldViscosity"; }

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
};

}

#endif
