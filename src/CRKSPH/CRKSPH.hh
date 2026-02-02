//---------------------------------Spheral++----------------------------------//
// CRKSPH -- The CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 21:52:29 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_CRKSPH_hh__
#define __Spheral_CRKSPH_hh__

#include "CRKSPH/CRKSPHBase.hh"
#include "Geometry/CellFaceFlag.hh"
#include "RK/RKCorrectionParams.hh"

#include <memory>

namespace Spheral {
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosityHandle;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
template<typename Dimension, typename Value, size_t numElements> class PairwiseField;
class FileIO;
}

namespace Spheral {

template<typename Dimension>
class CRKSPH: public CRKSPHBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using ThirdRankTensor = typename Dimension::ThirdRankTensor;
  using FourthRankTensor = typename Dimension::FourthRankTensor;
  using FifthRankTensor = typename Dimension::FifthRankTensor;
  using SymTensor = typename Dimension::SymTensor;
  using FacetedVolume = typename Dimension::FacetedVolume;

  using PairAccelerationsType = PairwiseField<Dimension, Vector, 1u>;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  CRKSPH(DataBase<Dimension>& dataBase,
                  ArtificialViscosityHandle<Dimension>& Q,
                  const RKOrder order,
                  const double cfl,
                  const bool useVelocityMagnitudeForDt,
                  const bool compatibleEnergyEvolution,
                  const bool evolveTotalEnergy,
                  const bool XSPH,
                  const MassDensityType densityUpdate,
                  const double epsTensile,
                  const double nTensile);

  // No default constructor, copying, or assignment.
  CRKSPH() = delete;
  CRKSPH(const CRKSPH&) = delete;
  CRKSPH& operator=(const CRKSPH&) = delete;

  // Destructor.
  virtual ~CRKSPH() = default;

  // Register the state Hydro expects to use and evolve.
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
                               const QType& Q) const;
    
  // The state field lists we're maintaining.
  const PairAccelerationsType& pairAccelerations() const { VERIFY2(mPairAccelerationsPtr, "CRKSPH ERROR: pairAccelerations not initialized on access"); return *mPairAccelerationsPtr; }

private:
  //--------------------------- Private Interface ---------------------------//
  std::unique_ptr<PairAccelerationsType> mPairAccelerationsPtr;
};

}

#endif
