//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
// References:
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//
// Created by JMO, Sun May 21 23:46:02 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_MonaghanGingoldViscosity__
#define __Spheral_MonaghanGingoldViscosity__

#include "ArtificialViscosity.hh"
#include "MonaghanGingoldViscosityView.hh"

namespace Spheral {

template<typename Dimension>
class MonaghanGingoldViscosity: public ArtificialViscosity<Dimension>,
                                public MonaghanGingoldViscosityView<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ViewType = MonaghanGingoldViscosityView<Dimension>;

  // Constructors.
  MonaghanGingoldViscosity(const Scalar Clinear,
                           const Scalar Cquadratic,
                           const TableKernel<Dimension>& kernel,
                           const bool linearInExpansion,
                           const bool quadraticInExpansion) :
    ArtificialViscosity<Dimension>(kernel),
    ViewType(Clinear,
             Cquadratic,
             linearInExpansion,
             quadraticInExpansion) {
  }

  virtual ~MonaghanGingoldViscosity() { }

  // No default construction, copying, or assignment
  MonaghanGingoldViscosity() = delete;
  MonaghanGingoldViscosity(const MonaghanGingoldViscosity&) = delete;
  MonaghanGingoldViscosity& operator=(const MonaghanGingoldViscosity&) = delete;

  // Restart methods.
  virtual std::string label()    const override { return "MonaghanGingoldViscosity"; }

  // Forward ArtificialViscosity's virtual parameter interface to ViewType
  // Each top level ArtificialViscosity value class must set these
  virtual Scalar Cl() const override { return ViewType::Cl(); }
  virtual Scalar Cq() const override { return ViewType::Cq(); }

  virtual void Cl(const Scalar x) override { ViewType::Cl(x); }
  virtual void Cq(const Scalar x) override { ViewType::Cq(x); }

  virtual bool balsaraShearCorrection() const override {
    return ViewType::balsaraShearCorrection();
  }
  virtual void balsaraShearCorrection(const bool x) override {
    ViewType::balsaraShearCorrection(x);
  }

  virtual Scalar epsilon2() const override {
    return ViewType::epsilon2();
  }
  virtual void epsilon2(const Scalar x) override {
    ViewType::epsilon2(x);
  }

  virtual Scalar negligibleSoundSpeed() const override {
    return ViewType::negligibleSoundSpeed();
  }
  virtual void negligibleSoundSpeed(const Scalar x) override {
    ViewType::negligibleSoundSpeed(x);
  }

  bool linearInExpansion() const {
    return ViewType::linearInExpansion();
  }
  bool quadraticInExpansion() const {
    return ViewType::quadraticInExpansion();
  }

  void linearInExpansion(const bool x) {
    ViewType::linearInExpansion(x);
  }
  void quadraticInExpansion(const bool x) {
    ViewType::quadraticInExpansion(x);
  }

  // View method
  ViewType view() const {
    return static_cast<const ViewType&>(*this);
  }
};

}

#endif
