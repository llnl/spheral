//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// Created by JMO, Thu Nov 20 14:13:18 PST 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_LimitedMonaghanGingoldViscosity__
#define __Spheral_LimitedMonaghanGingoldViscosity__

#include "LimitedMonaghanGingoldViscosityView.hh"

namespace Spheral {

template<typename Dimension>
class LimitedMonaghanGingoldViscosity : public MonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ArtViscView = ArtificialViscosityView<Dimension, Scalar>;
  using ViewType = LimitedMonaghanGingoldViscosityView<Dimension>;

  // Constructors.
  LimitedMonaghanGingoldViscosity(const Scalar Clinear,
                                  const Scalar Cquadratic,
                                  const TableKernel<Dimension>& kernel,
                                  const bool linearInExpansion,
                                  const bool quadraticInExpansion,
                                  const Scalar etaCritFrac,
                                  const Scalar etaFoldFrac) :
    MonaghanGingoldViscosity<Dimension>(Clinear, Cquadratic, kernel, 
                                        linearInExpansion, quadraticInExpansion),
    mEtaCritFrac(etaCritFrac),
    mEtaFoldFrac(etaFoldFrac) {
    m_viewPtr = chai::make_managed<ViewType>(Clinear,
                                             Cquadratic,
                                             linearInExpansion,
                                             quadraticInExpansion,
                                             etaCritFrac,
                                             etaFoldFrac);
  }
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
