//---------------------------------Spheral++----------------------------------//
// A simple form for the artificial viscosity due to Monaghan & Gingold.
// References:
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//
// Created by JMO, Sun May 21 23:46:02 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_MonaghanGingoldViscosityView__
#define __Spheral_MonaghanGingoldViscosityView__

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

  // Data access
  SPHERAL_HOST_DEVICE
  bool linearInExpansion() const { return mLinearInExpansion; }
  SPHERAL_HOST_DEVICE
  bool quadraticInExpansion() const { return mQuadraticInExpansion; }

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtificialViscosity (particularly DvDx).
  SPHERAL_HOST_DEVICE
  virtual void QPiij(Scalar& QPiij, Scalar& QPiji,      // result for QPi (Q/rho^2)
                     Scalar& Qij, Scalar& Qji,          // result for viscous pressure
                     const size_t nodeListi, const size_t i,
                     const size_t nodeListj, const size_t j,
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
protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mLinearInExpansion = false;
  bool mQuadraticInExpansion = false;

  using ArtificialViscosityBase<Dimension>::mClinear;
  using ArtificialViscosityBase<Dimension>::mCquadratic;
  using ArtificialViscosityBase<Dimension>::mEpsilon2;
  using ArtificialViscosityBase<Dimension>::mBalsaraShearCorrection;
};

}

#include "MonaghanGingoldViscosityViewInline.hh"

#endif
