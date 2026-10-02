//---------------------------------Spheral++----------------------------------//
// SolidSPH -- The SPH/ASPH solid material hydrodynamic package for Spheral++.
//
// Created by JMO, Fri Jul 30 11:07:33 PDT 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidSPH_RAJA_hh__
#define __Spheral_SolidSPH_RAJA_hh__

#include <float.h>
#include <string>

#include "SPH/SolidSPH.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
template<typename Dimension, typename Value, size_t numElements> class PairwiseField;
class FileIO;

template<typename Dimension>
class SolidSPH_RAJA: public SolidSPH<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using PairAccelerationsType = PairwiseField<Dimension, Vector, 1u>;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SolidSPH_RAJA(DataBase<Dimension>& dataBase,
                ArtificialViscosity<Dimension>& Q,
                const TableKernel<Dimension>& W,
                const TableKernel<Dimension>& WPi,
                const TableKernel<Dimension>& WGrad,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool useNewAccelerationMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool evolveTotalEnergy,
                const bool gradhCorrection,
                const bool XSPH,
                const bool correctVelocityGradient,
                const bool sumMassDensityOverAllNodeLists,
                const MassDensityType densityUpdate,
                const double epsTensile,
                const double nTensile,
                const bool damageRelieveRubble,
                const bool strengthInDamage,
                const Vector& xmin,
                const Vector& xmax);

  // No default constructor, copying, or assignment.
  SolidSPH_RAJA() = delete;
  SolidSPH_RAJA(const SolidSPH_RAJA&) = delete;
  SolidSPH_RAJA& operator=(const SolidSPH_RAJA&) = delete;

  // Destructor.
  virtual ~SolidSPH_RAJA() = default;

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
                               chai::managed_ptr<QType>& Q) const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  using SolidSPH<Dimension>::mDamageRelieveRubble;
  using SolidSPH<Dimension>::mStrengthInDamage;
};

}

#endif
