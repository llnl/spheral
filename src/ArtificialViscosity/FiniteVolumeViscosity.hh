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
class FiniteVolumeViscosity: public ArtificialViscosity<Dimension, typename Dimension::Scalar> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructor, destructor
  FiniteVolumeViscosity(const Scalar Clinear,
                        const Scalar Cquadratic,
                        const TableKernel<Dimension>& WT);
  virtual ~FiniteVolumeViscosity() = default;

  // No default construction, copying, or assignment
  FiniteVolumeViscosity() = delete;
  FiniteVolumeViscosity(const FiniteVolumeViscosity&) = delete;
  FiniteVolumeViscosity& operator=(const FiniteVolumeViscosity&) const = delete;

  // We are going to use a velocity gradient
  virtual bool requireVelocityGradient()             const override { return true; }

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

  // Override the method of computing the velocity gradient
  virtual void updateVelocityGradient(const DataBase<Dimension>& db,
                                      const State<Dimension>& state,
                                      const StateDerivatives<Dimension>& derivs) override;

  // Restart methods.
  virtual std::string label()                        const override { return "FiniteVolumeViscosity"; }

private:
  //--------------------------- Public Interface ---------------------------//
  using ArtificialViscosity<Dimension, Scalar>::mClinear;
  using ArtificialViscosity<Dimension, Scalar>::mCquadratic;
  using ArtificialViscosity<Dimension, Scalar>::mEpsilon2;
  using ArtificialViscosity<Dimension, Scalar>::mBalsaraShearCorrection;
};

}

#endif
