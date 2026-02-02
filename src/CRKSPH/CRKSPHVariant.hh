//---------------------------------Spheral++----------------------------------//
// CRKSPHVariant -- A development variant of CRKSPH for experimentation.
//
// Created by JMO, Thu Oct 12 14:24:43 PDT 2017
//----------------------------------------------------------------------------//
#ifndef __Spheral_CRKSPHVariant_hh__
#define __Spheral_CRKSPHVariant_hh__

#include <string>

#include "Physics/GenericHydro.hh"
#include "CRKSPHHydroBase.hh"

namespace Spheral {
  template<typename Dimension> class State;
  template<typename Dimension> class StateDerivatives;
  template<typename Dimension> class SmoothingScaleBase;
  template<typename Dimension> class ArtificialViscosity;
  template<typename Dimension> class TableKernel;
  template<typename Dimension> class DataBase;
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension, typename DataType> class FieldList;
  class FileIO;
}

namespace Spheral {

template<typename Dimension>
class CRKSPHVariant: public CRKSPHHydroBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  CRKSPHVariant(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                ArtificialViscosity<Dimension>& Q,
                const TableKernel<Dimension>& W,
                const TableKernel<Dimension>& WPi,
                const double filter,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool evolveTotalEnergy,
                const bool XSPH,
                const MassDensityType densityUpdate,
                const HEvolutionType HUpdate,
                const RKOrder correctionOrder,
                const RKVolumeType volumeType,
                const double epsTensile,
                const double nTensile,
                const bool limitMultimaterialTopology);

  // Destructor.
  virtual ~CRKSPHVariant() override;

  // An optional hook to initialize once when the problem is starting up.
  // Typically this is used to size arrays once all the materials and NodeLists have
  // been created.  It is assumed after this method has been called it is safe to
  // call Physics::registerState for instance to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;

  // Initialize the Hydro before we start a derivative evaluation.
  virtual
  bool initialize(const Scalar time,
                  const Scalar dt,
                  const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) override;
                          

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  CRKSPHVariant();
  CRKSPHVariant(const CRKSPHVariant&);
  CRKSPHVariant& operator=(const CRKSPHVariant&);
};

}

#endif
