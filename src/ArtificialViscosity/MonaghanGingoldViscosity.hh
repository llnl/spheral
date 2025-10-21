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

#include "ArtificialViscosityView.hh"
#include "ArtificialViscosity.hh"
#include "Utilities/CHAI_MA_wrapper.hh"

namespace Spheral {

// Forward declare so it can be made a friend of the view class
template<typename Dimension> class MonaghanGingoldViscosity;

template<typename Dimension>
class MonaghanGingoldViscosityView:
    public ArtificialViscosityView<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors.
  SPHERAL_HOST_DEVICE
  MonaghanGingoldViscosityView(const Scalar Clinear,
                               const Scalar Cquadratic,
                               const bool linearInExpansion,
                               const bool quadraticInExpansion) :
    ArtificialViscosityView<Dimension, Scalar>(Clinear,
                                               Cquadratic),
    mLinearInExpansion(linearInExpansion),
    mQuadraticInExpansion(quadraticInExpansion) {}

  SPHERAL_HOST_DEVICE virtual ~MonaghanGingoldViscosityView() = default;

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtificialViscosity (particularly DvDx).
  //SPHERAL_HOST_DEVICE
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
                     const FieldList<Dimension, Scalar>& fCl,
                     const FieldList<Dimension, Scalar>& fCq,
                     const FieldList<Dimension, Tensor>& DvDx) const override;

  friend class ArtificialViscosity<Dimension>;
  friend class MonaghanGingoldViscosity<Dimension>;
protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mLinearInExpansion = false;
  bool mQuadraticInExpansion = false;

  using ArtificialViscosityBase<Dimension>::mClinear;
  using ArtificialViscosityBase<Dimension>::mCquadratic;
  using ArtificialViscosityBase<Dimension>::mEpsilon2;
  using ArtificialViscosityBase<Dimension>::mBalsaraShearCorrection;
};

template<typename Dimension>
class MonaghanGingoldViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ArtViscView = ArtificialViscosityView<Dimension, Scalar>;
  using ViewType = MonaghanGingoldViscosityView<Dimension>;

  // Constructors.
  MonaghanGingoldViscosity(const Scalar Clinear,
                           const Scalar Cquadratic,
                           const TableKernel<Dimension>& kernel,
                           const bool linearInExpansion,
                           const bool quadraticInExpansion);

  virtual ~MonaghanGingoldViscosity() { m_viewPtr.free(); }

  // No default construction, copying, or assignment
  MonaghanGingoldViscosity() = delete;
  MonaghanGingoldViscosity(const MonaghanGingoldViscosity&) = delete;
  MonaghanGingoldViscosity& operator=(const MonaghanGingoldViscosity&) = delete;

  virtual std::type_index QPiTypeIndex() const override {
    return std::type_index(typeid(Scalar));
  }

  virtual chai::managed_ptr<ArtViscView> getScalarView() const override {
    return chai::dynamic_pointer_cast<ArtViscView, ViewType>(m_viewPtr);
  }

  // Restart methods.
  virtual std::string label()    const override { return "MonaghanGingoldViscosity"; }

  // Access data members
  bool linearInExpansion()                const { return mLinearInExpansion; }
  bool quadraticInExpansion()             const { return mQuadraticInExpansion; }
  void linearInExpansion(const bool x)          { mLinearInExpansion = x; updateManagedPtr(); }
  void quadraticInExpansion(const bool x)       { mQuadraticInExpansion = x; updateManagedPtr(); }

  // New member variables mLinearInExpansion and mQuadraticInExpansion require
  // this
  template<typename ViewPtr>
  void updateMembers(chai::managed_ptr<ViewPtr> a_viewPtr) {
    ArtificialViscosity<Dimension>::updateMembers(a_viewPtr);
    ASSIGN_MEMBER_ALL(m_viewPtr, mLinearInExpansion, mLinearInExpansion);
    ASSIGN_MEMBER_ALL(m_viewPtr, mQuadraticInExpansion, mQuadraticInExpansion);
  }

  virtual void updateManagedPtr() override {
    updateMembers(m_viewPtr);
  }
protected:
  // Not ideal but there is repeated member data between the value and view
  bool mLinearInExpansion = false;
  bool mQuadraticInExpansion = false;
private:
  std::type_index m_viewType = typeid(ViewType);
  chai::managed_ptr<ViewType> m_viewPtr;
};

}

#endif
