//---------------------------------Spheral++----------------------------------//
// A finite-volume based viscosity.  Assumes you have constructed the
// tessellation in the state.
//
// Created by JMO, Tue Aug 13 09:43:37 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_FiniteVolumeViscosity__
#define __Spheral_FiniteVolumeViscosity__

#include "ArtificialViscosity.hh"

namespace Spheral {

template<typename Dimension>
class FiniteVolumeViscosityView:
    public ArtificialViscosityView<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor
  SPHERAL_HOST_DEVICE
  FiniteVolumeViscosityView(const Scalar Clinear, const Scalar Cquadratic) :
    ArtificialViscosityView<Dimension, Scalar>(Clinear, Cquadratic) {}

  SPHERAL_HOST_DEVICE
  virtual ~FiniteVolumeViscosityView() = default;

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtficialViscosity (particularly DvDx).
  SPHERAL_HOST_DEVICE
  virtual void QPiij(Scalar& QPiij, Scalar& QPiji,      // result for QPi (Q/rho^2)
                     Scalar& Qij, Scalar& Qji,          // result for viscous pressure
                     const unsigned nodeListi, const unsigned i,
                     const unsigned nodeListj, const unsigned j,
                     const Vector& xi,
                     const SymTensor& Hi,
                     const Vector& etai,
                     const Vector& vi,
                     const Scalar rhoi,
                     const Scalar csi,
                     const Vector& xj,
                     const SymTensor& Hj,
                     const Vector& etaj,
                     const Vector& vj,
                     const Scalar rhoj,
                     const Scalar csj,
                     const FieldListView<Dimension, Scalar>& fCl,
                     const FieldListView<Dimension, Scalar>& fCq,
                     const FieldListView<Dimension, Tensor>& DvDx) const override;
  friend class ArtificialViscosity<Dimension>;
protected:
  //--------------------------- Protected Interface ---------------------------//
  using ArtificialViscosityBase<Dimension>::mClinear;
  using ArtificialViscosityBase<Dimension>::mCquadratic;
  using ArtificialViscosityBase<Dimension>::mEpsilon2;
  using ArtificialViscosityBase<Dimension>::mBalsaraShearCorrection;
  using ArtificialViscosityBase<Dimension>::mNegligibleSoundSpeed;
};

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
                        const TableKernel<Dimension>& WT);

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

  virtual chai::managed_ptr<ArtViscView> getScalarView() const override {
    return chai::dynamic_pointer_cast<ArtViscView, ViewType>(m_viewPtr);
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
