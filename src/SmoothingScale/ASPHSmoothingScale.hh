//---------------------------------Spheral++----------------------------------//
// ASPHSmoothingScale
//
// Implements the ASPH tensor smoothing scale algorithm.
//
// Created by JMO, Wed Sep 14 15:01:13 PDT 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_ASPHSmooothingScale__
#define __Spheral_ASPHSmooothingScale__

#include "SmoothingScale/SmoothingScaleBase.hh"
#include "Utilities/Functors.hh"

#include <utility>  // std::pair
#include <tuple>

namespace Spheral {

template<typename Dimension>
class ASPHSmoothingScale: public SmoothingScaleBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using FacetedVolume = typename Dimension::FacetedVolume;
  using HidealFilterType = typename SmoothingScaleBase<Dimension>::HidealFilterType;
  using RadialFunctorType = typename SmoothingScaleBase<Dimension>::RadialFunctorType;

  // Constructors, destructor.
  ASPHSmoothingScale(const HEvolutionType HUpdate,
                     const TableKernel<Dimension>& W,
                     const bool fixShape = false,
                     const bool radialOnly = false);
  ASPHSmoothingScale() = delete;
  virtual ~ASPHSmoothingScale() = default;

  // An optional hook to initialize once when the problem is starting up.
  // This is called after the materials and NodeLists are created. This method
  // should set the sizes of all arrays owned by the physics package and initialize
  // independent variables.
  // It is assumed after this method has been called it is safe to call
  // Physics::registerState to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

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

 // Similarly packages might want a hook to do some post-step finalizations.
  // Really we should rename this post-step finalize.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs) override;

  // Apply boundary conditions to the physics specific fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;

  // We require the Voronoi-like cells per point
  virtual bool requireVoronoiCells() const override                               { return (this->HEvolution() == HEvolutionType::IdealH and
                                                                                            not (mFixShape or mRadialOnly)); }

  // Access our internal data
  const TableKernel<Dimension>&                          WT()            const    { return mWT; }
  const FieldList<Dimension, Scalar>&                    zerothMoment()  const    { return mZerothMoment; }
  const FieldList<Dimension, SymTensor>&                 secondMoment()  const    { return mSecondMoment; }
  const FieldList<Dimension, SymTensor>&                 cellSecondMoment() const { return mCellSecondMoment; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "ASPHSmoothingScale"; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  const TableKernel<Dimension>& mWT;
  FieldList<Dimension, Scalar> mZerothMoment;
  FieldList<Dimension, SymTensor> mSecondMoment, mCellSecondMoment;

  using SmoothingScaleBase<Dimension>::mFixShape;
  using SmoothingScaleBase<Dimension>::mRadialOnly;
  using SmoothingScaleBase<Dimension>::mHidealFilterPtr;
  using SmoothingScaleBase<Dimension>::mRadialFunctorPtr;
  using SmoothingScaleBase<Dimension>::mRadius0;
};

}

#endif
