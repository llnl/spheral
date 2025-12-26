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
class TensorMonaghanGingoldViscosity : public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ArtViscView = ArtificialViscosityView<Dimension, Tensor>;
  using ViewType = TensorMonaghanGingoldViscosityView<Dimension>;

  // Constructors and destuctor
  TensorMonaghanGingoldViscosity(const Scalar Clinear,
                                 const Scalar Cquadratic,
                                 const TableKernel<Dimension>& kernel) :
    ArtificialViscosity<Dimension>(Clinear, Cquadratic, kernel) {
    m_viewPtr = chai::make_managed<TensorMonaghanGingoldViscosityView<Dimension>>(Clinear,
                                                                                  Cquadratic);
  }

  virtual ~TensorMonaghanGingoldViscosity() { m_viewPtr.free(); }

  // No default construction, copying, or assignment
  TensorMonaghanGingoldViscosity() = delete;
  TensorMonaghanGingoldViscosity(const TensorMonaghanGingoldViscosity&) = delete;
  TensorMonaghanGingoldViscosity& operator=(const TensorMonaghanGingoldViscosity&) = delete;

  // We need the velocity gradient
  virtual bool requireVelocityGradient() const override { return true; }

  // Restart methods.
  virtual std::string label() const override { return "TensorMonaghanGingoldViscosity"; }

  // View methods
  virtual std::type_index QPiTypeIndex() const override {
    return std::type_index(typeid(Tensor));
  }

  virtual chai::managed_ptr<ArtViscView> getTensorView() const override {
    return chai::dynamic_pointer_cast<ArtViscView, ViewType>(m_viewPtr);
  }

  // Useful for testing
  chai::managed_ptr<ViewType> getView() const {
    return m_viewPtr;
  }
protected:
  //--------------------------- Protected Interface ---------------------------//
  // Can simplify this call because no new value/view member data is made
  virtual void updateManagedPtr() override { this->updateMembers(m_viewPtr); }

private:
  std::type_index m_viewType = typeid(ViewType);
  chai::managed_ptr<ViewType> m_viewPtr;
};

}

#endif
