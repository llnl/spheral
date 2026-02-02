//---------------------------------Spheral++----------------------------------//
// SPHSmoothingScale
//
// Implements the standard SPH scalar smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 14:55:16 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_SPHSmooothingScale__
#define __Spheral_SPHSmooothingScale__

#include "SmoothingScale/SmoothingScaleBase.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

template<typename Dimension>
class SPHSmoothingScale: public SmoothingScaleBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor.
  SPHSmoothingScale(const HEvolutionType HUpdate,
                    const TableKernel<Dimension>& W);
  SPHSmoothingScale() = delete;
  virtual ~SPHSmoothingScale() {}

  // An optional hook to initialize once when the problem is starting up.
  // This is called after the materials and NodeLists are created. This method
  // should set the sizes of all arrays owned by the physics package and initialize
  // independent variables.
  // It is assumed after this method has been called it is safe to call
  // Physics::registerState to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Increment the derivatives.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  // Access our internal data
  const TableKernel<Dimension>&                          WT()            const { return mWT; }
  const FieldList<Dimension, Scalar>&                    zerothMoment()  const { return mZerothMoment; }
  const FieldList<Dimension, Vector>&                    firstMoment()   const { return mFirstMoment; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SPHSmoothingScale"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mWT;
  FieldList<Dimension, Scalar> mZerothMoment;
  FieldList<Dimension, Vector> mFirstMoment;
};

}

#endif
