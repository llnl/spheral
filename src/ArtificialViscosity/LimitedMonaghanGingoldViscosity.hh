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

template<typename Dimension>
class LimitedMonaghanGingoldViscosity: public MonaghanGingoldViscosity<Dimension> {
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
  LimitedMonaghanGingoldViscosity(const Scalar Clinear,
                                  const Scalar Cquadratic,
                                  const TableKernel<Dimension>& kernel,
                                  const bool linearInExpansion,
                                  const bool quadraticInExpansion,
                                  const Scalar etaCritFrac,
                                  const Scalar etaFoldFrac);
  virtual ~LimitedMonaghanGingoldViscosity() = default;

  // No default construction, copying, or assignment
  LimitedMonaghanGingoldViscosity() = delete;
  LimitedMonaghanGingoldViscosity(const LimitedMonaghanGingoldViscosity&) = delete;
  LimitedMonaghanGingoldViscosity& operator=(const LimitedMonaghanGingoldViscosity&) = delete;

  // We need the velocity gradient
  virtual bool requireVelocityGradient() const override { return true; }

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtficialViscosity (particularly DvDx).
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

  // Access our data
  Scalar etaCritFrac()                       const { return mEtaCritFrac; }
  Scalar etaFoldFrac()                       const { return mEtaFoldFrac; }

  void etaCritFrac(const Scalar x)                 { mEtaCritFrac = x; }
  void etaFoldFrac(const Scalar x)                 { mEtaFoldFrac = x; }

  // Restart methods.
  virtual std::string label()       const override { return "LimitedMonaghanGingoldViscosity"; }

protected:
  //--------------------------- Private Interface ---------------------------//
  double mEtaCritFrac, mEtaFoldFrac;

  using MonaghanGingoldViscosity<Dimension>::mLinearInExpansion;
  using MonaghanGingoldViscosity<Dimension>::mQuadraticInExpansion;
  using ArtificialViscosityView<Dimension, Scalar>::mClinear;
  using ArtificialViscosityView<Dimension, Scalar>::mCquadratic;
  using ArtificialViscosityView<Dimension, Scalar>::mEpsilon2;
  using ArtificialViscosityView<Dimension, Scalar>::mBalsaraShearCorrection;
};

}

#endif
