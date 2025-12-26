//---------------------------------Spheral++----------------------------------//
// HydroFieldNames -- A collection of standard Field names for the hydro 
// physics package.
//
// Created by JMO, Sat Aug 28 21:16:14 2004
//----------------------------------------------------------------------------//
#ifndef _Spheral_HydroFieldNames_
#define _Spheral_HydroFieldNames_

#include <string>

namespace Spheral {

struct HydroFieldNames {
  static inline const std::string mass = "mass";
  static inline const std::string position = "position";
  static inline const std::string velocity = "velocity";
  static inline const std::string H = "H";
  static inline const std::string work = "work";
  static inline const std::string velocityGradient = "velocity gradient";
  static inline const std::string internalVelocityGradient = "internal velocity gradient";
  // Non-hydro (things that don't modify the thermal energy) use this
  static inline const std::string acceleration = "delta " + velocity;   // Note here we *must* start with "delta " to work with IncrementFieldList!
  // Normal hydro sources (things that do modify material thermal energy)
  static inline const std::string hydroAcceleration = acceleration + " hydro";
  static inline const std::string ahgAcceleration = "delta " + hydroAcceleration + " anti hourglass";
  static inline const std::string massDensity = "mass density";
  static inline const std::string normalization = "normalization";
  static inline const std::string specificThermalEnergy = "specific thermal energy";
  static inline const std::string maxViscousPressure = "max viscous pressure";
  static inline const std::string effectiveViscousPressure = "effective viscous pressure";
  static inline const std::string massDensityCorrection = "density summation correction";
  static inline const std::string XSPHDeltaV = "XSPH delta vi";
  static inline const std::string XSPHWeightSum = "XSPH weight sum";
  static inline const std::string Hsmooth = "H smooth";
  static inline const std::string massZerothMoment = "mass zeroth moment";
  static inline const std::string massFirstMoment = "mass first moment";
  static inline const std::string massSecondMoment = "mass second moment";
  static inline const std::string pressure = "pressure";
  static inline const std::string partialPpartialEps = "partial pressure partial eps energy derivative";
  static inline const std::string partialPpartialRho = "partial pressure partial rho derivative";
  static inline const std::string temperature = "temperature";
  static inline const std::string soundSpeed = "sound speed";
  static inline const std::string pairAccelerations = "pair-wise accelerations";
  static inline const std::string pairWork = "pair-wise work";
  static inline const std::string selfAccelerations = "self-accelerations";
  static inline const std::string omegaGradh = "grad h corrections";
  static inline const std::string gamma = "ratio of specific heats";
  static inline const std::string entropy = "entropy";
  static inline const std::string PSPHcorrection = "PSPH Correction";
  static inline const std::string numberDensitySum = "number density sum";
  static inline const std::string timeStepMask = "time step mask";
  static inline const std::string surfacePoint = "surface point";
  static inline const std::string voidPoint = "void point";
  static inline const std::string etaVoidPoints = "eta void points";
  static inline const std::string cells = "cells";
  static inline const std::string cellFaceFlags = "cell face flags";
  static inline const std::string M_SPHCorrection = "M SPH gradient correction";
  static inline const std::string volume = "node volume";
  static inline const std::string linearMomentum = "linear momentum";
  static inline const std::string totalEnergy = "total energy";
  static inline const std::string mesh = "mesh";
  static inline const std::string hourglassMask = "hourglass mask";
  static inline const std::string faceVelocity = "face velocity";
  static inline const std::string faceForce = "face force";
  static inline const std::string faceMass = "face mass";
  static inline const std::string polyvols = "poly faceted volumes";
  static inline const std::string massDensityGradient = "mass density gradient";
  static inline const std::string ArtificialViscousClMultiplier = "Cl multiplier for artificial viscosity";
  static inline const std::string ArtificialViscousCqMultiplier = "Cq multiplier for artificial viscosity";
  static inline const std::string ArtificialViscosityVelocityGradient = velocityGradient + " for artificial viscosity";
  static inline const std::string pairQPi = "Pairwise artificial viscosity (P/rho^2)";
  static inline const std::string specificHeat = "specific heat";
  static inline const std::string normal = "outward normal direction";
  static inline const std::string surfaceArea = "boundary surface area";
};

}

#endif
