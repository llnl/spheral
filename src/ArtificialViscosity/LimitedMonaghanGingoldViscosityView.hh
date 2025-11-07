//---------------------------------Spheral++----------------------------------//
// Modified form of the standard SPH pair-wise viscosity due to Monaghan &
// Gingold.  This form is modified to use the velocity gradient to limit the
// velocity jump at the mid-point between points.
//
// Created by LDO, Thu Nov 6 13:23:00 PST 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_LimitedMonaghanGingoldViscosityView__
#define __Spheral_LimitedMonaghanGingoldViscosityView__

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

}

#endif
