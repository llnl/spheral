//---------------------------------Spheral++----------------------------------//
// ArtificialViscosityView -- The view class for all ArtificialViscosities in 
// Spheral++.
//
// This class contains the pure virtual function QPiij
// and specifies the QPi = Q/rho^2 return type,
// generally either a Scalar or a Tensor.
//
// Created by JMO, Sun May 21 21:16:43 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_ArtificialViscosityView__
#define __Spheral_ArtificialViscosityView__

#include "config.hh"
#include "chai/managed_ptr.hpp"
#include "Field/FieldList.hh"
#include "ArtificialViscosityBase.hh"

#include <utility>
#include <typeindex>

namespace Spheral {

template<typename Dimension, typename QPiType>
class ArtificialViscosityView : public ArtificialViscosityBase<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor
  SPHERAL_HOST_DEVICE
  ArtificialViscosityView(const Scalar Clinear,
                          const Scalar Cquadratic) :
    ArtificialViscosityBase<Dimension>(Clinear, Cquadratic) {}

  SPHERAL_HOST_DEVICE
  virtual ~ArtificialViscosityView() = default;

  // No default constructor, copying, or assignment
  // ArtificialViscosityView() = delete;
  // ArtificialViscosityView(const ArtificialViscosityView&) = delete;
  // ArtificialViscosityView& operator=(const ArtificialViscosityView&) = delete;

  //...........................................................................
  // Virtual methods we expect ArtificialViscosities to provide

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtficialViscosity (particularly DvDx).
  //SPHERAL_HOST_DEVICE
  virtual void QPiij(QPiType& QPiij, QPiType& QPiji,    // result for QPi (Q/rho^2)
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
                     const FieldList<Dimension, Tensor>& DvDx) const = 0;
protected:
  //--------------------------- Protected Interface ---------------------------//
  using ArtificialViscosityBase<Dimension>::mClinear;
  using ArtificialViscosityBase<Dimension>::mCquadratic;
  using ArtificialViscosityBase<Dimension>::mEpsilon2;
  using ArtificialViscosityBase<Dimension>::mBalsaraShearCorrection;
  using ArtificialViscosityBase<Dimension>::mNegligibleSoundSpeed;
};

}

#endif
