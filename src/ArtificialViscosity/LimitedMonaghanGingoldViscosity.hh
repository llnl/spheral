//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_LimitedMonaghanGingoldViscosity__
#define __Spheral_LimitedMonaghanGingoldViscosity__

#include "MonaghanGingoldViscosity.hh"

namespace Spheral {

// Forward declare so it can be made a friend of the view class
template<typename Dimension> class LimitedMonaghanGingoldViscosity;

template<typename Dimension>
class LimitedMonaghanGingoldViscosityView
  : public MonaghanGingoldViscosityView<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;
  using FourthRankTensor = typename Dimension::FourthRankTensor;
  using FifthRankTensor = typename Dimension::FifthRankTensor;

  // Constructors.
  SPHERAL_HOST_DEVICE
  LimitedMonaghanGingoldViscosityView(const Scalar Clinear,
                                      const Scalar Cquadratic,
                                      const bool linearInExpansion,
                                      const bool quadraticInExpansion,
                                      const Scalar etaCritFrac,
                                      const Scalar etaFoldFrac) :
    MonaghanGingoldViscosityView<Dimension>(Clinear, Cquadratic,
                                            linearInExpansion, quadraticInExpansion),
    mEtaCritFrac(etaCritFrac),
    mEtaFoldFrac(etaFoldFrac) {}

  SPHERAL_HOST_DEVICE virtual ~LimitedMonaghanGingoldViscosityView() = default;

  // Data access
  SPHERAL_HOST_DEVICE
  Scalar etaCritFrac() const { return mEtaCritFrac; }
  SPHERAL_HOST_DEVICE
  Scalar etaFoldFrac() const { return mEtaFoldFrac; }

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtficialViscosity (particularly DvDx).
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
                     const FieldListView<Dimension, Scalar>& fCl,
                     const FieldListView<Dimension, Scalar>& fCq,
                     const FieldListView<Dimension, Tensor>& DvDx) const override;

  friend class ArtificialViscosity<Dimension>;
  friend class MonaghanGingoldViscosity<Dimension>;
  friend class LimitedMonaghanGingoldViscosity<Dimension>;
protected:
  //--------------------------- Protected Interface ---------------------------//
  Scalar mEtaCritFrac;
  Scalar mEtaFoldFrac;

  using MonaghanGingoldViscosityView<Dimension>::mLinearInExpansion;
  using MonaghanGingoldViscosityView<Dimension>::mQuadraticInExpansion;
  using ArtificialViscosityBase<Dimension>::mClinear;
  using ArtificialViscosityBase<Dimension>::mCquadratic;
  using ArtificialViscosityBase<Dimension>::mEpsilon2;
  using ArtificialViscosityBase<Dimension>::mBalsaraShearCorrection;
};

template<typename Dimension>
class LimitedMonaghanGingoldViscosity : public MonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;
  using FourthRankTensor = typename Dimension::FourthRankTensor;
  using FifthRankTensor = typename Dimension::FifthRankTensor;
  using ArtViscView = ArtificialViscosityView<Dimension, Scalar>;
  using ViewType = LimitedMonaghanGingoldViscosityView<Dimension>;

  // Constructors.
  LimitedMonaghanGingoldViscosity(const Scalar Clinear,
                                  const Scalar Cquadratic,
                                  const TableKernel<Dimension>& kernel,
                                  const bool linearInExpansion,
                                  const bool quadraticInExpansion,
                                  const Scalar etaCritFrac,
                                  const Scalar etaFoldFrac);
  virtual ~LimitedMonaghanGingoldViscosity() { m_viewPtr.free(); }

  // No default construction, copying, or assignment
  LimitedMonaghanGingoldViscosity() = delete;
  LimitedMonaghanGingoldViscosity(const LimitedMonaghanGingoldViscosity&) = delete;
  LimitedMonaghanGingoldViscosity& operator=(const LimitedMonaghanGingoldViscosity&) = delete;

  // We need the velocity gradient
  virtual bool requireVelocityGradient() const override { return true; }

  // Access our data
  Scalar etaCritFrac()                 const { return mEtaCritFrac; }
  Scalar etaFoldFrac()                 const { return mEtaFoldFrac; }

  void etaCritFrac(const Scalar x)           { mEtaCritFrac = x; updateManagedPtr(); }
  void etaFoldFrac(const Scalar x)           { mEtaFoldFrac = x; updateManagedPtr(); }

  // Restart methods.
  virtual std::string label() const override { return "LimitedMonaghanGingoldViscosity"; }

  // View methods
  virtual std::type_index QPiTypeIndex() const override {
    return std::type_index(typeid(Scalar));
  }

  virtual chai::managed_ptr<ArtViscView> getScalarView() const override {
    return chai::dynamic_pointer_cast<ArtViscView, ViewType>(m_viewPtr);
  }

  // Useful for testing
  chai::managed_ptr<ViewType> getView() const {
    return m_viewPtr;
  }
protected:
  //--------------------------- Protected Interface ---------------------------//
  template<typename ViewPtr>
  void updateMembers(chai::managed_ptr<ViewPtr> a_viewPtr) {
    MonaghanGingoldViscosity<Dimension>::updateMembers(a_viewPtr);
    ASSIGN_MEMBER_ALL(a_viewPtr, mEtaCritFrac, mEtaCritFrac);
    ASSIGN_MEMBER_ALL(a_viewPtr, mEtaFoldFrac, mEtaFoldFrac);
  }
    
  virtual void updateManagedPtr() override { updateMembers(m_viewPtr); }

  // Not ideal but there is repeated member data between the value and view
  Scalar mEtaCritFrac = 1.0;
  Scalar mEtaFoldFrac = 0.2;
  using MonaghanGingoldViscosity<Dimension>::mLinearInExpansion;
  using MonaghanGingoldViscosity<Dimension>::mQuadraticInExpansion;
private:
  std::type_index m_viewType = typeid(ViewType);
  chai::managed_ptr<ViewType> m_viewPtr;
};

}

#endif
