//---------------------------------Spheral++----------------------------------//
// A finite-volume based viscosity.  Assumes you have constructed the
// tessellation in the state.
//
// Created by JMO, Tue Aug 13 09:43:37 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_FiniteVolumeViscosity__
#define __Spheral_FiniteVolumeViscosity__

#include "FiniteVolumeViscosityView.hh"

namespace Spheral {

template<typename Dimension>
class FiniteVolumeViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ArtViscView = ArtificialViscosityView<Dimension, Scalar>;
  using ViewType = FiniteVolumeViscosityView<Dimension>;

  // Constructor, destructor
  FiniteVolumeViscosity(const Scalar Clinear,
                        const Scalar Cquadratic,
                        const TableKernel<Dimension>& WT) :
    ArtificialViscosity<Dimension>(Clinear, Cquadratic, WT) { }

  virtual ~FiniteVolumeViscosity() { m_viewPtr.free(); }

  // No default construction, copying, or assignment
  FiniteVolumeViscosity() = delete;
  FiniteVolumeViscosity(const FiniteVolumeViscosity&) = delete;
  FiniteVolumeViscosity& operator=(const FiniteVolumeViscosity&) const = delete;

  // We are going to use a velocity gradient
  virtual bool requireVelocityGradient()             const override { return true; }

  // Override the method of computing the velocity gradient
  virtual void updateVelocityGradient(const DataBase<Dimension>& db,
                                      const State<Dimension>& state,
                                      const StateDerivatives<Dimension>& derivs) override;

  // Restart methods.
  virtual std::string label()                        const override { return "FiniteVolumeViscosity"; }

  // View methods
  virtual std::type_index QPiTypeIndex() const override {
    return std::type_index(typeid(Scalar));
  }

  virtual chai::managed_ptr<ArtViscView> getScalarView() override {
    return chai::dynamic_pointer_cast<ArtViscView, ViewType>(m_viewPtr);
  }

  // Useful for testing
  chai::managed_ptr<ViewType> getView() {
    initView();
    return m_viewPtr;
  }

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Initialize the managed pointer if it doesn't exist
  void initView() {
    if (!m_viewPtr) {
      m_viewPtr = chai::make_managed<ViewType>(mClinear, mCquadratic);
    }
  }

  // Reinitialize the managed pointer if it exists so member variables are up to date
  virtual void updateManagedPtr() override {
    if (m_viewPtr) {
      m_viewPtr.free();
      initView();
    }
  }
  using ArtificialViscosity<Dimension>::mClinear;
  using ArtificialViscosity<Dimension>::mCquadratic;
private:
  std::type_index m_viewType = typeid(ViewType);
  chai::managed_ptr<ViewType> m_viewPtr = nullptr;
};

}

#endif
