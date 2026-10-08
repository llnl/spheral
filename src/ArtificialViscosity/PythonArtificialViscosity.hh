//---------------------------------Spheral++----------------------------------//
// PythonArtificialViscosity -- A Python-overridable artificial viscosity
// for rapid prototyping of new viscosity models.
//
// This class allows users to implement custom artificial viscosity algorithms
// in Python without C++ compilation. It uses a simplified interface with 12
// parameters instead of the full 20+ parameter QPiij signature.
//
// WARNING: CPU-only, significantly slower than C++ implementations (~100-1000x).
// Use for prototyping only, not production runs.
//
// Created by JMO and some agents, September 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_PythonArtificialViscosity__
#define __Spheral_PythonArtificialViscosity__

#include "ArtificialViscosity.hh"
#include "ArtificialViscosityView.hh"

#include <tuple>
#include <type_traits>

namespace Spheral {

template<typename Dimension, typename QPiType>
class PythonArtificialViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ArtViscViewScalar = ArtificialViscosityView<Dimension, Scalar>;
  using ArtViscViewTensor = ArtificialViscosityView<Dimension, Tensor>;
  using ArtViscView = ArtificialViscosityView<Dimension, QPiType>;

  // Constructor
  PythonArtificialViscosity(const Scalar Clinear,
                            const Scalar Cquadratic,
                            const TableKernel<Dimension>& kernel);

  virtual ~PythonArtificialViscosity();

  // No default constructor, copying, or assignment
  PythonArtificialViscosity() = delete;
  PythonArtificialViscosity(const PythonArtificialViscosity&) = delete;
  PythonArtificialViscosity& operator=(const PythonArtificialViscosity&) = delete;

  //...........................................................................
  // SIMPLIFIED virtual method for Python to override (CPU-only)
  virtual std::tuple<QPiType, QPiType, Scalar, Scalar>                                   // Outputs: (Qij/rho^2, Qji/rho^2, Qij, Qji)
  computeQPiij(const Vector& xi, const Vector& vi, const Scalar rhoi, const Scalar csi,  // Particle i
               const Vector& xj, const Vector& vj, const Scalar rhoj, const Scalar csj,  // Particle j
               const Vector& etai, const Vector& etaj,                                   // Pre-computed H*(xi-xj) and H*(xj-xi)
               const Scalar fCli, const Scalar fCqi,                                     // Pre-computed multipliers for i
               const Scalar fClj, const Scalar fCqj) const = 0;                          // Pre-computed multipliers for j

  //...........................................................................
  // Standard ArtificialViscosity interface

  // Return the appropriately typed CPU-only View wrapper.
  virtual chai::managed_ptr<ArtViscViewScalar> getScalarView() override;
  virtual chai::managed_ptr<ArtViscViewTensor> getTensorView() override;

  // Return the QPi output type.
  virtual std::type_index QPiTypeIndex() const override {
    return std::type_index(typeid(QPiType));
  }

  // Label for restart
  virtual std::string label() const override { return "PythonArtificialViscosity"; }

protected:
  //--------------------------- Protected Interface ---------------------------//
  // No managed ptr to update (View wrapper is different pattern)
  virtual void updateManagedPtr() override {}

private:
  //--------------------------- Private Interface ---------------------------//
  // Forward declare the private View class
  class PythonAVView;

  // Managed pointer to View wrapper (created on-demand)
  chai::managed_ptr<PythonAVView> mView;
};

//------------------------------------------------------------------------------
// Private View implementation that wraps Python AV
// This class implements the full QPiij signature and calls back to the
// simplified computeQPiij method on the Python class
//------------------------------------------------------------------------------
template<typename Dimension, typename QPiType>
class PythonArtificialViscosity<Dimension, QPiType>::PythonAVView
    : public ArtificialViscosityView<Dimension, QPiType> {
public:
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructor - store pointer to parent Python AV
  PythonAVView(PythonArtificialViscosity<Dimension, QPiType>* parent);

  virtual ~PythonAVView() = default;

  //...........................................................................
  // Implement full QPiij signature by extracting from FieldListView
  // and calling simplified computeQPiij
  virtual void QPiij(QPiType& QPiij, QPiType& QPiji,
                     Scalar& Qij, Scalar& Qji,
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

private:
  // Pointer back to parent Python AV class
  PythonArtificialViscosity<Dimension, QPiType>* mParent;
};

}

#endif
