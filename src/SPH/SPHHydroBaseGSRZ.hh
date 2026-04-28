//---------------------------------Spheral++----------------------------------//
// SPHHydroBaseGSRZ -- The SPH/ASPH hydrodynamic package for Spheral++,
//                     specialized for 2D RZ (cylindrical) geometry.
//
// Based on the methodology described in the dissertation
// de Catalunya Departament de Física i Enginyeria Nuclear, U. P., García Senz, D., (2012).
// AxisSPH:devising and validating an axisymmetric smoothed particle hydrodynamics code.
// TDX (Tesis Doctorals en Xarxa). Universitat Politècnica de Catalunya.
//
// Note this version is currently abusing our ordinary 2D geometric types,
// implicitly mapping x->z, y->r.
//
// Created by JMO, Tue Apr 26 16:06:21 PDT 2016
//----------------------------------------------------------------------------//
#ifndef __Spheral_SPHHydroBaseGSRZ_hh__
#define __Spheral_SPHHydroBaseGSRZ_hh__

#include "SPHHydroBase.hh"
#include "Geometry/Dimension.hh"

#include <string>

namespace Spheral {

class SPHHydroBaseGSRZ: public SPHHydroBase<Dim<2> > {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef Dim<2> Dimension;
  typedef Dimension::Scalar Scalar;
  typedef Dimension::Vector Vector;
  typedef Dimension::Tensor Tensor;
  typedef Dimension::SymTensor SymTensor;

  typedef Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  // Constructors.
  SPHHydroBaseGSRZ(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                   DataBase<Dimension>& dataBase,
                   ArtificialViscosityView<Dimension>& Q,
                   const TableKernel<Dimension>& W,
                   const TableKernel<Dimension>& WPi,
                   const double filter,
                   const double cfl,
                   const bool useVelocityMagnitudeForDt,
                   const bool compatibleEnergyEvolution,
                   const bool evolveTotalEnergy,
                   const bool gradhCorrection,
                   const bool XSPH,
                   const bool correctVelocityGradient,
                   const bool sumMassDensityOverAllNodeLists,
                   const MassDensityType densityUpdate,
                   const HEvolutionType HUpdate,
                   const double epsTensile,
                   const double nTensile,
                   const Vector& xmin,
                   const Vector& xmax);

  // Destructor.
  virtual ~SPHHydroBaseGSRZ();

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Evaluate the derivatives for the principle hydro variables:
  // mass density, velocity, and specific thermal energy.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;

  // This method is called once at the beginning of a timestep, after all state registration.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;
               
  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SPHHydroBaseGSRZ" ; }
  //****************************************************************************

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  SPHHydroBaseGSRZ();
  SPHHydroBaseGSRZ(const SPHHydroBaseGSRZ&);
  SPHHydroBaseGSRZ& operator=(const SPHHydroBaseGSRZ&);
};

}

#endif
