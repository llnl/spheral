//---------------------------------Spheral++----------------------------------//
// HydroFieldNames -- A collection of standard Field names for the hydro 
// physics package.
//
// Created by JMO, Sat Aug 28 21:16:14 2004
//----------------------------------------------------------------------------//

#include "HydroFieldNames.hh"

const std::string Spheral::HydroFieldNames::mass = "mass";
const std::string Spheral::HydroFieldNames::position = "position";
const std::string Spheral::HydroFieldNames::velocity = "velocity";
const std::string Spheral::HydroFieldNames::H = "H";
const std::string Spheral::HydroFieldNames::work = "work";
const std::string Spheral::HydroFieldNames::velocityGradient = "velocity gradient";
const std::string Spheral::HydroFieldNames::internalVelocityGradient = "internal velocity gradient";
const std::string Spheral::HydroFieldNames::acceleration = "delta " + Spheral::HydroFieldNames::velocity;   // Note here we *must* start with "delta " to work with IncrementFieldList!
const std::string Spheral::HydroFieldNames::hydroAcceleration = Spheral::HydroFieldNames::acceleration + " hydro";
const std::string Spheral::HydroFieldNames::ahgAcceleration = "delta " + Spheral::HydroFieldNames::hydroAcceleration + " anti hourglass";
const std::string Spheral::HydroFieldNames::massDensity = "mass density";
const std::string Spheral::HydroFieldNames::normalization = "normalization";
const std::string Spheral::HydroFieldNames::specificThermalEnergy = "specific thermal energy";
const std::string Spheral::HydroFieldNames::maxViscousPressure = "max viscous pressure";
const std::string Spheral::HydroFieldNames::effectiveViscousPressure = "effective viscous pressure";
const std::string Spheral::HydroFieldNames::massDensityCorrection = "density summation correction";
const std::string Spheral::HydroFieldNames::XSPHDeltaV = "XSPH delta vi";
const std::string Spheral::HydroFieldNames::XSPHWeightSum = "XSPH weight sum";
const std::string Spheral::HydroFieldNames::Hsmooth = "H smooth";
const std::string Spheral::HydroFieldNames::massZerothMoment = "mass zeroth moment";
const std::string Spheral::HydroFieldNames::massFirstMoment = "mass first moment";
const std::string Spheral::HydroFieldNames::massSecondMoment = "mass second moment";
const std::string Spheral::HydroFieldNames::pressure = "pressure";
const std::string Spheral::HydroFieldNames::partialPpartialEps = "partial pressure partial eps energy derivative";
const std::string Spheral::HydroFieldNames::partialPpartialRho = "partial pressure partial rho derivative";
const std::string Spheral::HydroFieldNames::temperature = "temperature";
const std::string Spheral::HydroFieldNames::soundSpeed = "sound speed";
const std::string Spheral::HydroFieldNames::pairAccelerations = "pair-wise accelerations";
const std::string Spheral::HydroFieldNames::pairWork = "pair-wise work";
const std::string Spheral::HydroFieldNames::selfAccelerations = "self-accelerations";
const std::string Spheral::HydroFieldNames::omegaGradh = "grad h corrections";
const std::string Spheral::HydroFieldNames::gamma = "ratio of specific heats";
const std::string Spheral::HydroFieldNames::entropy = "entropy";
const std::string Spheral::HydroFieldNames::PSPHcorrection = "PSPH Correction";
const std::string Spheral::HydroFieldNames::numberDensitySum = "number density sum";
const std::string Spheral::HydroFieldNames::timeStepMask = "time step mask";
const std::string Spheral::HydroFieldNames::surfacePoint = "surface point";
const std::string Spheral::HydroFieldNames::voidPoint = "void point";
const std::string Spheral::HydroFieldNames::etaVoidPoints = "eta void points";
const std::string Spheral::HydroFieldNames::cells = "cells";
const std::string Spheral::HydroFieldNames::cellFaceFlags = "cell face flags";
const std::string Spheral::HydroFieldNames::M_SPHCorrection = "M SPH gradient correction";
const std::string Spheral::HydroFieldNames::volume = "node volume";
const std::string Spheral::HydroFieldNames::linearMomentum = "linear momentum";
const std::string Spheral::HydroFieldNames::totalEnergy = "total energy";
const std::string Spheral::HydroFieldNames::mesh = "mesh";
const std::string Spheral::HydroFieldNames::hourglassMask = "hourglass mask";
const std::string Spheral::HydroFieldNames::faceVelocity = "face velocity";
const std::string Spheral::HydroFieldNames::faceForce = "face force";
const std::string Spheral::HydroFieldNames::faceMass = "face mass";
const std::string Spheral::HydroFieldNames::polyvols = "poly faceted volumes";
const std::string Spheral::HydroFieldNames::massDensityGradient = "mass density gradient";
const std::string Spheral::HydroFieldNames::ArtificialViscousClMultiplier = "Cl multiplier for artificial viscosity";
const std::string Spheral::HydroFieldNames::ArtificialViscousCqMultiplier = "Cq multiplier for artificial viscosity";
const std::string Spheral::HydroFieldNames::ArtificialViscosityVelocityGradient = Spheral::HydroFieldNames::velocityGradient + " for artificial viscosity";
const std::string Spheral::HydroFieldNames::pairQPi = "Pairwise artificial viscosity (P/rho^2)";
const std::string Spheral::HydroFieldNames::specificHeat = "specific heat";
const std::string Spheral::HydroFieldNames::normal = "outward normal direction";
const std::string Spheral::HydroFieldNames::surfaceArea = "boundary surface area";
