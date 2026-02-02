//---------------------------------Spheral++----------------------------------//
// SphericalSPH -- An SPH/ASPH hydrodynamic package for Spheral++ specialized
//                 for 1D Spherical (r) geometry.
//
// Based on the algorithm described in
// Omang, M., Børve, S., & Trulsen, J. (2006). SPH in spherical and cylindrical coordinates.
// Journal of Computational Physics, 213(1), 391412. https://doi.org/10.1016/j.jcp.2005.08.023
//
// Note this version is currently abusing our ordinary 1D geometric types,
// implicitly mapping x->r.
//
// Created by JMO, Tue Dec 22 10:04:21 PST 2020
//----------------------------------------------------------------------------//
#ifndef __Spheral_SphericalSPH_hh__
#define __Spheral_SphericalSPH_hh__

#include "SPH/SPHBase.hh"
#include "Kernel/SphericalKernel.hh"
#include "Geometry/Dimension.hh"

#include <utility>   // pair
#include <string>
#include <memory>    // unique_ptr

namespace Spheral {

class SphericalSPH: public SPHBase<Dim<1>> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<1>;
  using Scalar = Dimension::Scalar;
  using Vector = Dimension::Vector;
  using Tensor = Dimension::Tensor;
  using SymTensor = Dimension::SymTensor;

  using PairAccelerationsType = PairwiseField<Dimension, Vector, 2u>;
  using ConstBoundaryIterator = Physics<Dimension>::ConstBoundaryIterator;

  // Constructors.
  SphericalSPH(DataBase<Dimension>& dataBase,
               ArtificialViscosityHandle<Dimension>& Q,
               const SphericalKernel& W,
               const SphericalKernel& WPi,
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
  SphericalSPH() = delete;
  SphericalSPH(const SphericalSPH&) = delete;
  SphericalSPH& operator=(const SphericalSPH&) = delete;

  // Destructor.
  virtual ~SphericalSPH() = default;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  // This method is called once at the beginning of a timestep, after all state registration.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
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
  template<typename QType>
  void evaluateDerivativesImpl(const Scalar time,
                               const Scalar dt,
                               const DataBase<Dimension>& dataBase,
                               const State<Dimension>& state,
                               StateDerivatives<Dimension>& derivatives,
                               const QType& Q) const;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;
               
  // Access the stored interpolation kernels.
  // These hide the base class "kernel" methods which return vanilla TableKernels.
  const SphericalKernel& kernel()                         const { return mKernel; }
  const SphericalKernel& PiKernel()                       const { return mPiKernel; }

  // We also have a funny self-Q term for interactions near the origin.
  double Qself()                                          const { return mQself; }
  void Qself(const double x)                                    { mQself = x; }

  // Access our state.
  const PairAccelerationsType& pairAccelerations()        const { VERIFY2(mPairAccelerationsPtr, "SPH ERROR: pairAccelerations not initialized on access"); return *mPairAccelerationsPtr; }
  const FieldList<Dimension, Vector>& selfAccelerations() const { return mSelfAccelerations; }

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label()                    const override { return "SphericalSPH" ; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  double mQself;

  // The specialized kernels
  const SphericalKernel& mKernel;
  const SphericalKernel& mPiKernel;

  std::unique_ptr<PairAccelerationsType> mPairAccelerationsPtr;
  FieldList<Dimension, Vector> mSelfAccelerations;
};

}

#endif
