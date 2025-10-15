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

namespace Spheral {

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
  //SPHERAL_HOST_DEVICE
  MonaghanGingoldViscosityView(const Scalar Clinear,
                               const Scalar Cquadratic,
                               const bool linearInExpansion,
                               const bool quadraticInExpansion) :
    ArtificialViscosityView<Dimension, Tensor>(Clinear, Cquadratic),
    mLinearInExpansion(linearInExpansion),
    mQuadraticInExpansion(quadraticInExpansion) {}

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


protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mLinearInExpansion = false;
  bool mQuadraticInExpansion = false;

  using ArtificialViscosityView<Dimension, Scalar>::mClinear;
  using ArtificialViscosityView<Dimension, Scalar>::mCquadratic;
  using ArtificialViscosityView<Dimension, Scalar>::mEpsilon2;
  using ArtificialViscosityView<Dimension, Scalar>::mBalsaraShearCorrection;
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
    return m_viewPtr;
  }

  // Restart methods.
  virtual std::string label()    const override { return "MonaghanGingoldViscosity"; }

  // Access data members
  bool linearInExpansion()                const { return m_viewPtr->mLinearInExpansion; }
  bool quadraticInExpansion()             const { return m_viewPtr->mQuadraticInExpansion; }
  void linearInExpansion(const bool x)          { m_viewPtr->mLinearInExpansion = x; }
  void quadraticInExpansion(const bool x)       { m_viewPtr->mQuadraticInExpansion = x; }
protected:
  std::type_index m_viewType = typeid(MonaghanGingoldViscosityView<Dimension>);
  chai::managed_ptr<ArtViscView> m_viewPtr;
};

}

#endif
