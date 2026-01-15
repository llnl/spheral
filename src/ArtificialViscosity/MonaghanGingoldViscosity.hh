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
#include "MonaghanGingoldViscosityView.hh"

namespace Spheral {

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
                           const bool quadraticInExpansion) :
    ArtificialViscosity<Dimension>(Clinear, Cquadratic, kernel),
    mLinearInExpansion(linearInExpansion),
    mQuadraticInExpansion(quadraticInExpansion) {
    m_viewPtr = chai::make_managed<MonaghanGingoldViscosityView<Dimension>>(Clinear,
                                                                            Cquadratic,
                                                                            linearInExpansion,
                                                                            quadraticInExpansion);
  }

  virtual ~MonaghanGingoldViscosity() { m_viewPtr.free(); }

  // No default construction, copying, or assignment
  MonaghanGingoldViscosity() = delete;
  MonaghanGingoldViscosity(const MonaghanGingoldViscosity&) = delete;
  MonaghanGingoldViscosity& operator=(const MonaghanGingoldViscosity&) = delete;

  // Restart methods.
  virtual std::string label()    const override { return "MonaghanGingoldViscosity"; }

  // Access data members
  bool linearInExpansion()                const { return mLinearInExpansion; }
  bool quadraticInExpansion()             const { return mQuadraticInExpansion; }
  void linearInExpansion(const bool x)          { mLinearInExpansion = x; updateManagedPtr(); }
  void quadraticInExpansion(const bool x)       { mQuadraticInExpansion = x; updateManagedPtr(); }

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
  // New member variables like mLinearInExpansion require this
  template<typename ViewPtr>
  void updateMembers(chai::managed_ptr<ViewPtr> a_viewPtr) {
    ArtificialViscosity<Dimension>::updateMembers(a_viewPtr);
    ASSIGN_MEMBER_ALL(a_viewPtr, mLinearInExpansion, mLinearInExpansion);
    ASSIGN_MEMBER_ALL(a_viewPtr, mQuadraticInExpansion, mQuadraticInExpansion);
  }

  virtual void updateManagedPtr() override { updateMembers(m_viewPtr); }

  // Not ideal but there is repeated member data between the value and view
  bool mLinearInExpansion = false;
  bool mQuadraticInExpansion = false;
private:
  std::type_index m_viewType = typeid(ViewType);
  chai::managed_ptr<ViewType> m_viewPtr;
};

}

#endif
