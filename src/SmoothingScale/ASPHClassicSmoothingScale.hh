//---------------------------------Spheral++----------------------------------//
// ASPHClassicSmoothingScale
//
// Implements our classic ASPH algorithm (2010 SPHERIC proceedings style)
//
// Created by JMO, Fri Jan  3 10:48:43 PST 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_ASPHClassicSmooothingScale__
#define __Spheral_ASPHClassicSmooothingScale__

#include "SmoothingScale/SmoothingScaleBase.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

template<typename Dimension>
class ASPHClassicSmoothingScale: public SmoothingScaleBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using HidealFilterType = typename SmoothingScaleBase<Dimension>::HidealFilterType;
  using RadialFunctorType = typename SmoothingScaleBase<Dimension>::RadialFunctorType;

  // Constructors, destructor.
  ASPHClassicSmoothingScale(const HEvolutionType HUpdate,
                            const TableKernel<Dimension>& W,
                            const bool fixShape = false,
                            const bool radialOnly = false);
  ASPHClassicSmoothingScale() = delete;
  virtual ~ASPHClassicSmoothingScale() = default;

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
  const FieldList<Dimension, SymTensor>&                 secondMoment()  const { return mSecondMoment; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "ASPHClassicSmoothingScale"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mWT;
  FieldList<Dimension, Scalar> mZerothMoment;
  FieldList<Dimension, Vector> mFirstMoment;
  FieldList<Dimension, SymTensor> mSecondMoment;

  using SmoothingScaleBase<Dimension>::mFixShape;
  using SmoothingScaleBase<Dimension>::mRadialOnly;
  using SmoothingScaleBase<Dimension>::mHidealFilterPtr;
  using SmoothingScaleBase<Dimension>::mRadialFunctorPtr;
  using SmoothingScaleBase<Dimension>::mRadius0;
};

}

#endif
