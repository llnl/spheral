//---------------------------------Spheral++----------------------------------//
// PythonScalarArtificialViscosity -- A Python-overridable artificial viscosity
// for rapid prototyping of new viscosity models.
//
// This class allows users to implement custom artificial viscosity algorithms
// in Python without C++ compilation. It uses a simplified interface with 12
// parameters instead of the full 20+ parameter QPiij signature.
//
// WARNING: CPU-only, significantly slower than C++ implementations (~100-1000x).
// Use for prototyping only, not production runs.
//
// Created by Claude, September 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_PythonArtificialViscosity__
#define __Spheral_PythonArtificialViscosity__

#include "ArtificialViscosity.hh"
#include "ArtificialViscosityView.hh"

namespace Spheral {

template<typename Dimension>
class PythonScalarArtificialViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ArtViscViewScalar = ArtificialViscosityView<Dimension, Scalar>;

  // Constructor
  PythonScalarArtificialViscosity(const Scalar Clinear,
                                  const Scalar Cquadratic,
                                  const TableKernel<Dimension>& kernel);

  virtual ~PythonScalarArtificialViscosity();

  // No default constructor, copying, or assignment
  PythonScalarArtificialViscosity() = delete;
  PythonScalarArtificialViscosity(const PythonScalarArtificialViscosity&) = delete;
  PythonScalarArtificialViscosity& operator=(const PythonScalarArtificialViscosity&) = delete;

  //...........................................................................
  // SIMPLIFIED virtual method for Python to override (CPU-only)
  // This method has only 12 parameters instead of the full 20+ QPiij signature
  // All FieldListView lookups and complex data are pre-computed by the C++ wrapper
  virtual void computeQPiij(Scalar& QPiij, Scalar& QPiji,                                             // Outputs: Q/rho^2
                            Scalar& Qij, Scalar& Qji,                                                 // Outputs: viscous pressure Q
                            const Vector& xi, const Vector& vi, const Scalar rhoi, const Scalar csi,  // Particle i
                            const Vector& xj, const Vector& vj, const Scalar rhoj, const Scalar csj,  // Particle j
                            const Vector& etai, const Vector& etaj,                                   // Pre-computed H*(xi-xj) and H*(xj-xi)
                            const Scalar fCli, const Scalar fCqi,                                     // Pre-computed multipliers for i
                            const Scalar fClj, const Scalar fCqj) const = 0;                          // Pre-computed multipliers for j

  //...........................................................................
  // Standard ArtificialViscosity interface

  // Override getScalarView to provide CPU-only View wrapper
  virtual chai::managed_ptr<ArtViscViewScalar> getScalarView() override;

  // Return proper type index (Scalar)
  virtual std::type_index QPiTypeIndex() const override {
    return std::type_index(typeid(Scalar));
  }

  // Label for restart
  virtual std::string label() const override { return "PythonScalarArtificialViscosity"; }

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
template<typename Dimension>
class PythonScalarArtificialViscosity<Dimension>::PythonAVView
    : public ArtificialViscosityView<Dimension, typename Dimension::Scalar> {
public:
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructor - store pointer to parent Python AV
  PythonAVView(PythonScalarArtificialViscosity<Dimension>* parent);

  virtual ~PythonAVView() = default;

  //...........................................................................
  // Implement full QPiij signature by extracting from FieldListView
  // and calling simplified computeQPiij
  virtual void QPiij(Scalar& QPiij, Scalar& QPiji,
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
  PythonScalarArtificialViscosity<Dimension>* mParent;
};

}

#endif
