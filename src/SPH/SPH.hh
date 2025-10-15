//---------------------------------Spheral++----------------------------------//
// SPH -- The classic SPH/ASPH hydrodynamic packages for Spheral++.
// 
// Created by JMO, Thu Nov 21 16:36:40 PST 2024
//----------------------------------------------------------------------------//
#ifndef __Spheral_SPH__
#define __Spheral_SPH__

#include "SPH/SPHBase.hh"
#include "Utilities/CHAI_MA_wrapper.hh"

#include <memory>

namespace Spheral {

template<typename Dimension> class Physics;
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
template<typename Dimension, typename Value, size_t numElements> class PairwiseField;

template<typename Dimension>
class SPH: public SPHBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using DimensionType = Dimension;
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using PairAccelerationsType = PairwiseField<Dimension, Vector, 1u>;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SPH(DataBase<Dimension>& dataBase,
      ArtificialViscosity<Dimension>& Q,
      const TableKernel<Dimension>& W,
      const TableKernel<Dimension>& WPi,
      const double cfl,
      const bool useVelocityMagnitudeForDt,
      const bool compatibleEnergyEvolution,
      const bool evolveTotalEnergy,
      const bool gradhCorrection,
      const bool XSPH,
      const bool correctVelocityGradient,
      const bool sumMassDensityOverAllNodeLists,
      const MassDensityType densityUpdate,
      const double epsTensile,
      const double nTensile,
      const Vector& xmin,
      const Vector& xmax);

  // No default constructor, copying, or assignment.
  SPH() = delete;
  SPH(const SPH&) = delete;
  SPH& operator=(const SPH&) = delete;

  // Destructor.
  virtual ~SPH() = default;

  // Register the state
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;
  template<typename QType>
  void evaluateDerivativesImpl(const Scalar time,
                               const Scalar dt,
                               const DataBase<Dimension>& dataBase,
                               const State<Dimension>& state,
                               StateDerivatives<Dimension>& derivatives,
                               chai::managed_ptr<QType> Q) const;

  // Access our state.
  const PairAccelerationsType& pairAccelerations() const { VERIFY2(mPairAccelerationsPtr, "SPH ERROR: pairAccelerations not initialized on access"); return *mPairAccelerationsPtr; }

private:
  //---------------------------  Private Interface ---------------------------//
  std::unique_ptr<PairAccelerationsType> mPairAccelerationsPtr;
};

}

#endif
