//---------------------------------Spheral++----------------------------------//
// A specialized form of the TensorMonaghanGingoldViscosity for use with CRKSPH.
//
// Created by J. Michael Owen, Wed Nov  5 23:51:31 PST 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorCRKSPHViscosity__
#define __Spheral_TensorCRKSPHViscosity__

#include "TensorMonaghanGingoldViscosity.hh"
#include "RK/RKCorrectionParams.hh"

namespace Spheral {

template<typename Dimension>
class TensorCRKSPHViscosity: public TensorMonaghanGingoldViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;

  // Constructors, destructor
  TensorCRKSPHViscosity(const Scalar Clinear,
                        const Scalar Cquadratic,
                        const TableKernel<Dimension>& kernel,
                        const RKOrder order = RKOrder::LinearOrder);
  virtual ~TensorCRKSPHViscosity() = default;

  // No default construction, copying, or assignment
  TensorCRKSPHViscosity();
  TensorCRKSPHViscosity(const TensorCRKSPHViscosity&);
  TensorCRKSPHViscosity& operator=(const TensorCRKSPHViscosity&) const;

  // We do require the velocity gradient
  virtual bool requireVelocityGradient()                const override { return true; }

  // Override the method for computing the velocity gradient
  virtual void updateVelocityGradient(const DataBase<Dimension>& db,
                                      const State<Dimension>& state,
                                      const StateDerivatives<Dimension>& derivs) override;

  // Override the method listing our RK requirements
  virtual std::set<RKOrder> requireReproducingKernels() const override { return std::set<RKOrder>({mOrder}); }

  RKOrder order()                                       const          { return mOrder; }
  void order(const RKOrder x)                                          { mOrder = x; }

  // Restart methods.
  virtual std::string label()                           const override { return "TensorCRKSPHViscosity"; }

private:
  //--------------------------- Private Interface ---------------------------//
  RKOrder mOrder;
};

}

#endif
