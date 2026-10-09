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
#include "ArtificialViscosityVariant.hh"
#include "Field/FieldList.hh"

#include <tuple>

namespace Spheral {

template<typename Dimension, typename QPiType> class PythonArtificialViscosity;

//------------------------------------------------------------------------------
// CPU-only adapter that invokes the Python-overridable viscosity method.
// It intentionally owns no viscosity parameters: those belong to its parent.
//------------------------------------------------------------------------------
template<typename Dimension, typename QPiType>
class PythonArtificialViscosityCallView {
public:
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ReturnType = QPiType;

  PythonArtificialViscosityCallView(
    const PythonArtificialViscosity<Dimension, QPiType>* parent);

  void QPiij(QPiType& QPiij, QPiType& QPiji,
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
             const FieldListView<Dimension, Tensor>& DvDx) const;

private:
  const PythonArtificialViscosity<Dimension, QPiType>* mParent;
};

template<typename Dimension, typename QPiType>
class PythonArtificialViscosity: public ArtificialViscosity<Dimension>,
                                 public ArtificialViscosityView<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ViewType = PythonArtificialViscosityCallView<Dimension, QPiType>;

  // Constructor
  PythonArtificialViscosity(const Scalar Clinear,
                            const Scalar Cquadratic,
                            const TableKernel<Dimension>& kernel);

  virtual ~PythonArtificialViscosity() = default;

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

  // Forward ArtificialViscosity's virtual parameter interface to the view.
  virtual Scalar Cl() const override { return ArtificialViscosityView<Dimension>::Cl(); }
  virtual Scalar Cq() const override { return ArtificialViscosityView<Dimension>::Cq(); }
  virtual bool balsaraShearCorrection() const override {
    return ArtificialViscosityView<Dimension>::balsaraShearCorrection();
  }
  virtual Scalar epsilon2() const override { return ArtificialViscosityView<Dimension>::epsilon2(); }
  virtual Scalar negligibleSoundSpeed() const override {
    return ArtificialViscosityView<Dimension>::negligibleSoundSpeed();
  }

  virtual void Cl(const Scalar x) override { ArtificialViscosityView<Dimension>::Cl(x); }
  virtual void Cq(const Scalar x) override { ArtificialViscosityView<Dimension>::Cq(x); }
  virtual void balsaraShearCorrection(const bool x) override {
    ArtificialViscosityView<Dimension>::balsaraShearCorrection(x);
  }
  virtual void epsilon2(const Scalar x) override { ArtificialViscosityView<Dimension>::epsilon2(x); }
  virtual void negligibleSoundSpeed(const Scalar x) override {
    ArtificialViscosityView<Dimension>::negligibleSoundSpeed(x);
  }

  // Return a CPU-only adapter that refers to this Python object.
  ViewType view() const { return ViewType(this); }
  ArtificialViscosityVariant<Dimension> variantType() const override {
    return this->view();
  }

  // Label for restart.
  virtual std::string label() const override { return "PythonArtificialViscosity"; }
};

}

#endif
