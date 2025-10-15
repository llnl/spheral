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
#include <utility>
#include <typeindex>

namespace Spheral {

template<typename Dimension> class ArtificialViscosity;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension, typename QPiType>
class ArtificialViscosityView {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using ReturnType = QPiType;

  // Constructors, destructor
  //SPHERAL_HOST_DEVICE
  ArtificialViscosityView(const Scalar Clinear,
                          const Scalar Cquadratic,
                          const bool   BalsaraShearCorrection = false,
                          const Scalar Epsilon2 = 1.0E-2,
                          const Scalar NegligibleSoundSpeed = 1.0E-10) :
    mClinear(Clinear),
    mCquadratic(Cquadratic),
    mBalsaraShearCorrection(BalsaraShearCorrection),
    mEpsilon2(Epsilon2),
    mNegligibleSoundSpeed(NegligibleSoundSpeed) {}
  //SPHERAL_HOST_DEVICE
  virtual ~ArtificialViscosityView() = default;

  // No default constructor, copying, or assignment
  ArtificialViscosityView() = delete;
  ArtificialViscosityView(const ArtificialViscosityView&) = delete;
  ArtificialViscosityView& operator=(const ArtificialViscosityView&) = delete;

  //...........................................................................
  // Virtual methods we expect ArtificialViscosities to provide
  // Required method returning the type_index of the descendant QPiType
  //std::type_index QPiTypeIndex() const { return std::type_index(typeid(QPiType)); }

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
  //...........................................................................
  // Methods
  // Calculate the curl of the velocity given the stress tensor.
  //SPHERAL_HOST_DEVICE
  Scalar curlVelocityMagnitude(const Tensor& DvDx) const;

  // Find the Balsara shear correction multiplier
  //SPHERAL_HOST_DEVICE
  Scalar calcBalsaraShearCorrection(const Tensor& DvDx,
                                    const SymTensor& H,
                                    const Scalar& cs) const;
  friend class ArtificialViscosity<Dimension>;
protected:
  Scalar mClinear;
  Scalar mCquadratic;
  // Switch for the Balsara shear correction.
  bool mBalsaraShearCorrection;
  // Parameters for the Q limiter.
  Scalar mEpsilon2;
  Scalar mNegligibleSoundSpeed;
};

}

#include "ArtificialViscosityViewInline.hh"

#endif
