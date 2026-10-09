//---------------------------------Spheral++----------------------------------//
// View class for the MonaghanGingoldViscosity.
// References:
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//
// Created by JMO, Sun May 21 23:46:02 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_MonaghanGingoldViscosityView__
#define __Spheral_MonaghanGingoldViscosityView__

#include "ArtificialViscosityView.hh"

namespace Spheral {

template<typename Dimension>
class MonaghanGingoldViscosityView:
    public ArtificialViscosityView<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ReturnType = Scalar;
  using ArtificialViscosityView<Dimension>::Cl;
  using ArtificialViscosityView<Dimension>::Cq;
  using ArtificialViscosityView<Dimension>::balsaraShearCorrection;
  using ArtificialViscosityView<Dimension>::epsilon2;
  using ArtificialViscosityView<Dimension>::negligibleSoundSpeed;

  // Constructors.
  SPHERAL_HOST_DEVICE
  MonaghanGingoldViscosityView(const Scalar Clinear,
                               const Scalar Cquadratic,
                               const bool linearInExpansion,
                               const bool quadraticInExpansion) :
    ArtificialViscosityView<Dimension>(Clinear,
                                       Cquadratic),
    mLinearInExpansion(linearInExpansion),
    mQuadraticInExpansion(quadraticInExpansion) {}

  SPHERAL_HOST_DEVICE ~MonaghanGingoldViscosityView() = default;

  // Data access
  SPHERAL_HOST_DEVICE
  bool linearInExpansion() const { return mLinearInExpansion; }
  SPHERAL_HOST_DEVICE
  bool quadraticInExpansion() const { return mQuadraticInExpansion; }
  SPHERAL_HOST_DEVICE
  void linearInExpansion(const bool x) { mLinearInExpansion = x; }
  SPHERAL_HOST_DEVICE
  void quadraticInExpansion(const bool x) { mQuadraticInExpansion = x; }

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtificialViscosity (particularly DvDx).
  SPHERAL_HOST_DEVICE
  void QPiij(Scalar& QPiij, Scalar& QPiji,      // result for QPi (Q/rho^2)
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
             const FieldListView<Dimension, Tensor>& DvDx) const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mLinearInExpansion = false;
  bool mQuadraticInExpansion = false;

  using ArtificialViscosityView<Dimension>::mClinear;
  using ArtificialViscosityView<Dimension>::mCquadratic;
  using ArtificialViscosityView<Dimension>::mBalsaraShearCorrection;
  using ArtificialViscosityView<Dimension>::mEpsilon2;
};

}

#include "MonaghanGingoldViscosityViewInline.hh"

#endif
