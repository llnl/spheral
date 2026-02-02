//---------------------------------Spheral++----------------------------------//
// ArtificialViscosity -- The base class for all ArtificialViscosities in 
// Spheral++.
//
// This intermediate class between ArtificialViscosityHandle and the actual
// artficial viscosity implementations specifies the QPi = Q/rho^2 return type,
// generally either a Scalar or a Tensor.
//
// Created by JMO, Sun May 21 21:16:43 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_ArtificialViscosity__
#define __Spheral_ArtificialViscosity__

#include "ArtificialViscosity/ArtificialViscosityHandle.hh"

#include <utility>
#include <typeindex>

namespace Spheral {

template<typename Dimension, typename QPiType>
class ArtificialViscosity: public ArtificialViscosityHandle<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using ReturnType = QPiType;

  // Constructors, destructor
  ArtificialViscosity(const Scalar Clinear,
                      const Scalar Cquadratic,
                      const TableKernel<Dimension>& kernel): ArtificialViscosityHandle<Dimension>(Clinear, Cquadratic, kernel) {}
  virtual ~ArtificialViscosity() = default;

  // No default constructor, copying, or assignment
  ArtificialViscosity() = delete;
  ArtificialViscosity(const ArtificialViscosity&) = delete;
  ArtificialViscosity& operator=(const ArtificialViscosity&) = delete;

  //...........................................................................
  // Virtual methods we expect ArtificialViscosities to provide
  // Required method returning the type_index of the descendant QPiType
  virtual std::type_index QPiTypeIndex() const override { return std::type_index(typeid(QPiType)); }

  // All ArtificialViscosities must provide the pairwise QPi term (pressure/rho^2)
  // Returns the pair values QPiij and QPiji by reference as the first two arguments.
  // Note the final FieldLists (fCl, fCQ, DvDx) should be the special versions registered
  // by the ArtficialViscosity (particularly DvDx).
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
};

}

#endif
