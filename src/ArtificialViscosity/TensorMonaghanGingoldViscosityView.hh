//---------------------------------Spheral++----------------------------------//
// A modified form of the Monaghan & Gingold viscosity, extended to tensor
// formalism.
//
// Created by J. Michael Owen, Mon Sep  2 14:45:35 PDT 2002
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorMonaghanGingoldViscosityView__
#define __Spheral_TensorMonaghanGingoldViscosityView__

#include "ArtificialViscosityView.hh"
#include "ArtificialViscosity.hh"

namespace Spheral {

template<typename Dimension>
class TensorMonaghanGingoldViscosityView :
    public ArtificialViscosityView<Dimension, typename Dimension::Tensor> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  SPHERAL_HOST_DEVICE
  TensorMonaghanGingoldViscosityView(const Scalar Clinear, const Scalar Cquadratic) :
    ArtificialViscosityView<Dimension, Tensor>(Clinear, Cquadratic) {}

  SPHERAL_HOST_DEVICE
  virtual ~TensorMonaghanGingoldViscosityView() = default;

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtficialViscosity (particularly DvDx).
  SPHERAL_HOST_DEVICE
  virtual void QPiij(Tensor& QPiij, Tensor& QPiji,      // result for QPi (Q/rho^2)
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
protected:
  //--------------------------- Protected Interface ---------------------------//
  using ArtificialViscosityBase<Dimension>::mClinear;
  using ArtificialViscosityBase<Dimension>::mCquadratic;
  using ArtificialViscosityBase<Dimension>::mEpsilon2;
  using ArtificialViscosityBase<Dimension>::mBalsaraShearCorrection;
};

}

#include "TensorMonaghanGingoldViscosityViewInline.hh"

#endif
