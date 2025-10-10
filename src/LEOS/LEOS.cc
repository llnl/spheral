//---------------------------------Spheral++----------------------------------//
// LEOS
//
// Wrap up the Livermore Equation Of State (LEOS) package as an option in 
// Spheral.
//
// Created by JMO, Wed May  8 16:29:41 PDT 2013
//----------------------------------------------------------------------------//
#include "LEOS/LEOS.hh"
#include "LEOS/LEOS_bundle.hh"
#include "Field/Field.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/safeInv.hh"

#include <LEOS.h>

#include <algorithm>
#include <filesystem>
#include <sstream>

using std::min;
using std::max;
using std::abs;
using std::cout;
using std::cerr;
using std::endl;

using namespace LEOS;   // This namespace being the same as our class name causes confusion without this

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Generic method to lookup an EOS field (no derivatives)
//------------------------------------------------------------------------------
template<typename ContainerType>
inline
void LEOSLookup(ContainerType& result,
                const ContainerType& massDensity,
                const ContainerType& specificThermalEnergy,
                const typename LEOS_bundle::name_t fname,
                const typename LEOS_bundle::eosnum_t eosnum,
                const typename LEOS_bundle::name_t dbname,
                const double resultConv,
                const double rhoConv,
                const double epsConv) {

  const auto n = result.size();
  REQUIRE(massDensity.size() == n and specificThermalEnergy.size() == n);

  if (n > 0u) {

    // Terrible hack!  We convert the Fields in-place to LEOS units, which
    // means casting away constness.  Boo!
    double *rhoPtr = const_cast<double*>(&massDensity[0]);
    double *epsPtr = const_cast<double*>(&specificThermalEnergy[0]);

    // Convert to LEOS units
    if (rhoConv != 1.0 or epsConv != 1.0) {
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        rhoPtr[i] /= rhoConv;
        epsPtr[i] /= epsConv;
      }
    }

    // Get the LEOS functions
    auto fepsPtr = LEOS_bundle::functionPtr(LEOS_Et, eosnum, dbname);
    auto funcPtr = LEOS_bundle::functionPtr(fname, eosnum, dbname);

    LEOS_LookupOptions lopt;
    lopt.calculateFunction(true).calculateDFDX(false).calculateDFDY(false);

    // Do the LEOS lookup
    std::vector<double> T(n);
    double* dummy = nullptr;
    fepsPtr->inverseEval(n, rhoPtr, epsPtr, &T[0], dummy, dummy, &lopt);
    funcPtr->eval(n, rhoPtr, &T[0], &result[0], dummy, dummy, &lopt);

    // Convert back to our runtime units
    if (rhoConv != 1.0 or epsConv != 1.0 or resultConv != 1.0) {
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        result[i] *= resultConv;
        rhoPtr[i] *= rhoConv;
        epsPtr[i] *= epsConv;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Generic method to lookup an EOS field (with derivatives)
//------------------------------------------------------------------------------
template<typename ContainerType>
inline
void LEOSLookup(ContainerType& result,
                ContainerType& dfdeps,
                ContainerType& dfdrho,
                const ContainerType& massDensity,
                const ContainerType& specificThermalEnergy,
                const typename LEOS_bundle::name_t fname,
                const typename LEOS_bundle::eosnum_t eosnum,
                const typename LEOS_bundle::name_t dbname,
                const double resultConv,
                const double rhoConv,
                const double epsConv) {

  const auto n = result.size();
  REQUIRE(massDensity.size() == n and
          specificThermalEnergy.size() == n and
          dfdeps.size() == n and
          dfdrho.size() == n);

  if (n > 0u) {

    // Terrible hack!  We convert the Fields in-place to LEOS units, which
    // means casting away constness.  Boo!
    double *rhoPtr = const_cast<double*>(&massDensity[0]);
    double *epsPtr = const_cast<double*>(&specificThermalEnergy[0]);

    // Convert to LEOS units
    if (rhoConv != 1.0 or epsConv != 1.0) {
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        rhoPtr[i] /= rhoConv;
        epsPtr[i] /= epsConv;
      }
    }

    // Get the LEOS function
    auto fepsPtr = LEOS_bundle::functionPtr(LEOS_Et, eosnum, dbname);
    auto funcPtr = LEOS_bundle::functionPtr(fname, eosnum, dbname);

    LEOS_LookupOptions lopt;
    lopt.calculateFunction(true).calculateDFDX(true).calculateDFDY(true);

    // Do the LEOS lookup
    std::vector<double> T(n), dedt(n), dedr(n);
    fepsPtr->inverseEval(n, rhoPtr, epsPtr, &T[0], &dedr[0], &dedt[0], &lopt);
    funcPtr->eval(n, rhoPtr, &T[0], &result[0], &dfdrho[0], &dfdeps[0], &lopt);

    // Convert back to our runtime units.  If nothing else we need to convert
    // dfdeps from temperature to spec energy, so gotta do this loop.
    const auto dfdeConv = resultConv*safeInv(epsConv);
    const auto dfdrConv = resultConv*safeInv(rhoConv);
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      result[i] *= resultConv;
      dfdeps[i] *= dfdeConv*safeInv(dedt[i]);
      dfdrho[i] *= dfdrConv;
      rhoPtr[i] *= rhoConv;
      epsPtr[i] *= epsConv;
    }
  }
}

//------------------------------------------------------------------------------
// Generic method to lookup a single value (no derivatives)
//------------------------------------------------------------------------------
inline
double LEOSLookupValue(const double& massDensity,
                       const double& specificThermalEnergy,
                       const typename LEOS_bundle::name_t fname,
                       const typename LEOS_bundle::eosnum_t eosnum,
                       const typename LEOS_bundle::name_t dbname,
                       const double resultConv,
                       const double rhoConv,
                       const double epsConv) {

  // Convert to LEOS units
  const L8DOUBLE rho = massDensity/rhoConv;
  const L8DOUBLE eps = specificThermalEnergy/epsConv;

  // Get the LEOS functions
  auto fepsPtr = LEOS_bundle::functionPtr(LEOS_Et, eosnum, dbname);
  auto funcPtr = LEOS_bundle::functionPtr(fname, eosnum, dbname);

  LEOS_LookupOptions lopt;
  lopt.calculateFunction(true).calculateDFDX(false).calculateDFDY(false);

  // Do the LEOS lookup
  L8DOUBLE T, fval, dummy;
  fepsPtr->inverseEval(rho, eps, T, dummy, dummy, &lopt);
  funcPtr->eval(rho, T, fval, dummy, dummy, &lopt);

  // Return in our units
  return resultConv * fval;
}

//------------------------------------------------------------------------------
// Generic method to lookup a single value (with derivatives)
//------------------------------------------------------------------------------
inline
void LEOSLookupValue(double& f,
                     double& dfdeps,
                     double& dfdrho,
                     const double& massDensity,
                     const double& specificThermalEnergy,
                     const typename LEOS_bundle::name_t fname,
                     const typename LEOS_bundle::eosnum_t eosnum,
                     const typename LEOS_bundle::name_t dbname,
                     const double resultConv,
                     const double rhoConv,
                     const double epsConv) {

  // Convert to LEOS units
  const L8DOUBLE rho = massDensity/rhoConv;
  const L8DOUBLE eps = specificThermalEnergy/epsConv;

  // Get the LEOS functions
  auto fepsPtr = LEOS_bundle::functionPtr(LEOS_Et, eosnum, dbname);
  auto funcPtr = LEOS_bundle::functionPtr(fname, eosnum, dbname);

  LEOS_LookupOptions lopt;
  lopt.calculateFunction(true).calculateDFDX(true).calculateDFDY(true);

  // Do the LEOS lookup
  L8DOUBLE T, dedt, dedr;
  fepsPtr->inverseEval(rho, eps, T, dedr, dedt, &lopt);
  funcPtr->eval(rho, T, f, dfdrho, dfdeps, &lopt);   // Note dfdeps is actually dfdT coming from this call -- have to convert

  // Return in our units
  f *= resultConv;
  dfdeps *= resultConv*safeInv(dedt*epsConv);
  dfdrho *= resultConv*safeInv(rhoConv);
}

} // anonymous

//------------------------------------------------------------------------------
// Construct for the given material number.
//------------------------------------------------------------------------------
template<typename Dimension>
LEOS<Dimension>::
LEOS(const int materialNumber,
     const PhysicalConstants& constants,
     const double externalPressure,
     const double minimumPressure,
     const double maximumPressure,
     const double minimumPressureDamage,
     const MaterialPressureMinType minPressureType,
     const std::string dbname,
     const std::string leosFileFormat,
     const double atomicWeight):
  SolidEquationOfState<Dimension>(1.0,    // dummy out the reference density initially.
                                  0.0,
                                  std::numeric_limits<double>::max(),
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minimumPressureDamage,
                                  minPressureType,
                                  externalPressure),
  mMaterialNumber(materialNumber),
  mDatabaseName(dbname),
  mRhoConv(1.0),
  mTconv(1.0),
  mPconv(1.0),
  mEconv(1.0),
  mCVconv(1.0),
  mVelConv(1.0),
  mSconv(1.0),
  mAtomicWeight(atomicWeight),
  mK0(0.0),
  mLEOSunits(0.01,   // cm expressed as meters.
             0.001,  // g expressed in kg.
             1.0) {  // sec in secs.

  // Build our unit conversion factors.
  const double lconv = mLEOSunits.unitLengthMeters() / constants.unitLengthMeters(),
               mconv = mLEOSunits.unitMassKg() / constants.unitMassKg(),
               tconv = mLEOSunits.unitTimeSec() / constants.unitTimeSec();
  mRhoConv = mconv/(lconv*lconv*lconv);
  mPconv = mconv/(lconv*tconv*tconv);
  mTconv = 1.0; // Assuming LEOS reports temperatures in K
  mEconv = FastMath::square(lconv/tconv);
  mCVconv = mEconv/mTconv;
  mVelConv = lconv/tconv;
  mSconv = mEconv/(mconv*mTconv);

  // Check if we're using a special LEOS file format
  const std::unordered_map<std::string, fileType> LEOS_fileTypes = {{""            , fileType::PDBFILE_LEOS},
                                                                    {"PDB"         , fileType::PDBFILE_LEOS},
                                                                    {"HDF"         , fileType::HDF5FILE_LEOS},
                                                                    {"Sesame"      , fileType::SESAMEFILE_LEOS},
                                                                    {"OpacServer"  , fileType::OPACITY_SERVER_LEOS},
                                                                    {"Sesame*"     , fileType::SESAME_ASCII_LEOS},
                                                                    {"Hyades"      , fileType::HYADES_LEOS},
                                                                    {"Purgatorio"  , fileType::PURGATORIO_LEOS},
                                                                    {"ascii"       , fileType::ASCII_LEOS},
                                                                    {"OneD"        , fileType::ONEDTEXT_LEOS},
                                                                    {"Sample_LEOS" , fileType::SAMPLE_TEXT_LEOS}};
  const auto dbformat = LEOS_fileTypes.at(leosFileFormat);

  // Add this material and its functions to the LEOS data base
  VERIFY2(dbname == "leos" or
          std::filesystem::exists(dbname),
          "LEOS Error: unable to access data base file " << dbname);
  auto matPtr = LEOS_bundle::materialPtr(materialNumber, dbname, dbformat);

  // Set the actual reference density in the SolidEquationOfState.
  L8DOUBLE rho0;
  matPtr->getMetadata(MAT_RHO_0, rho0);
  SolidEquationOfState<Dimension>::referenceDensity(rho0*mRhoConv);

  // Set the atomic weight.
  if (mAtomicWeight <= 0.0) {
    std::vector<L8DOUBLE> speciesFrac, speciesWeight;
    matPtr->getMetadata(MAT_ATOMIC_FRACTION, speciesFrac);
    matPtr->getMetadata(MAT_ATOMIC_WEIGHT, speciesWeight);
    const auto nspecies = speciesFrac.size();
    VERIFY(nspecies > 0u);
    VERIFY(speciesWeight.size() == nspecies);
    VERIFY(fuzzyEqual(std::accumulate(speciesFrac.begin(), speciesFrac.end(), 0.0), 1.0, 1.0e-10));
    auto minFrac = 1.0;
    for (auto k = 0u; k < nspecies; ++k) {
      if (speciesFrac[k] > 0.0) minFrac = std::min(minFrac, speciesFrac[k]);
      mAtomicWeight += speciesFrac[k]*speciesWeight[k];
    }
    CHECK(minFrac > 0.0 and minFrac <= 1.0);
    mAtomicWeight /= minFrac;
  }
  CHECK(mAtomicWeight > 0.0);

  // And the reference bulk modulus
  matPtr->getMetadata(MAT_BULK_MODULUS, mK0);
  mK0 *= mPconv;

  ENSURE(this->valid());
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LEOS<Dimension>::
setPressure(Field<Dimension, Scalar>& pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  
  // Do the LEOS operation
  LEOSLookup(pressure, massDensity, specificThermalEnergy,
             LEOS_Pt, mMaterialNumber, mDatabaseName,
             mPconv, mRhoConv, mEconv);

  // Apply pressure limits
  const auto n = massDensity.size();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    pressure(i) = this->applyPressureLimits(pressure[i]);
  }
}

//------------------------------------------------------------------------------
// Set the pressure and derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
LEOS<Dimension>::
setPressureAndDerivs(Field<Dimension, Scalar>& pressure,
                     Field<Dimension, Scalar>& dPdu,
                     Field<Dimension, Scalar>& dPdrho,
                     const Field<Dimension, Scalar>& massDensity,
                     const Field<Dimension, Scalar>& specificThermalEnergy) const {

  const auto n = massDensity.size();
  if (n > 0) {

    // Do the LEOS operation
    LEOSLookup(pressure, dPdu, dPdrho,
               massDensity, specificThermalEnergy,
               LEOS_Pt, mMaterialNumber, mDatabaseName,
               mPconv, mRhoConv, mEconv);

#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      pressure(i) = this->applyPressureLimits(pressure(i));
    }
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LEOS<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {

  const auto n = temperature.size();
  if (n > 0u) {

    // Terrible hack!  We convert the Fields in-place to LEOS units, which
    // means casting away constness.  Boo!
    double *rhoPtr = const_cast<double*>(&massDensity[0]);
    double *epsPtr = const_cast<double*>(&specificThermalEnergy[0]);

    // Convert to LEOS units
    if (mRhoConv != 1.0 or mEconv != 1.0) {
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        rhoPtr[i] /= mRhoConv;
        epsPtr[i] /= mEconv;
      }
    }

    // Get the LEOS functions
    auto fepsPtr = LEOS_bundle::functionPtr(LEOS_Et, mMaterialNumber, mDatabaseName);

    LEOS_LookupOptions lopt;
    lopt.calculateFunction(true).calculateDFDX(false).calculateDFDY(false);

    // Do the LEOS lookup
    std::vector<double> dedt(1u), dedr(1u);
    fepsPtr->inverseEval(n, rhoPtr, epsPtr, &temperature[0], &dedr[0], &dedt[0], &lopt);

    // Convert back to our runtime units
    if (mRhoConv != 1.0 or mEconv != 1.0 or mTconv != 1.0) {
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        temperature[i] *= mTconv;
        rhoPtr[i] *= mRhoConv;
        epsPtr[i] *= mEconv;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LEOS<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {

  const auto n = specificThermalEnergy.size();
  if (n > 0u) {

    // Terrible hack!  We convert the Fields in-place to LEOS units, which
    // means casting away constness.  Boo!
    double *rhoPtr = const_cast<double*>(&massDensity[0]);
    double *Tptr = const_cast<double*>(&temperature[0]);

    // Convert to LEOS units
    if (mRhoConv != 1.0 or mTconv != 1.0) {
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        rhoPtr[i] /= mRhoConv;
        Tptr[i] /= mTconv;
      }
    }

    // Get the LEOS functions
    auto fepsPtr = LEOS_bundle::functionPtr(LEOS_Et, mMaterialNumber, mDatabaseName);

    LEOS_LookupOptions lopt;
    lopt.calculateFunction(true).calculateDFDX(false).calculateDFDY(false);

    // Do the LEOS lookup
    double* dummy = nullptr;
    fepsPtr->eval(n, rhoPtr, Tptr, &specificThermalEnergy[0], dummy, dummy, &lopt);

    // Convert back to our runtime units
    if (mRhoConv != 1.0 or mEconv != 1.0 or mTconv != 1.0) {
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        specificThermalEnergy[i] *= mEconv;
        rhoPtr[i] *= mRhoConv;
        Tptr[i] *= mTconv;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LEOS<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {

  const auto n = temperature.size();
  if (n > 0u) {

    // Terrible hack!  We convert the Fields in-place to LEOS units, which
    // means casting away constness.  Boo!
    double *rhoPtr = const_cast<double*>(&massDensity[0]);
    double *Tptr = const_cast<double*>(&temperature[0]);

    // Convert to LEOS units
    if (mRhoConv != 1.0 or mTconv != 1.0) {
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        rhoPtr[i] /= mRhoConv;
        Tptr[i] /= mTconv;
      }
    }

    // Get the LEOS functions
    auto fepsPtr = LEOS_bundle::functionPtr(LEOS_Et, mMaterialNumber, mDatabaseName);

    LEOS_LookupOptions lopt;
    lopt.calculateFunction(true).calculateDFDX(true).calculateDFDY(true);

    // Do the LEOS lookup
    std::vector<double> eps(n), dedr(n);
    fepsPtr->eval(n, rhoPtr, Tptr, &eps[0], &specificHeat[0], &dedr[0], &lopt);

    // Convert back to our runtime units
    if (mRhoConv != 1.0 or mTconv != 1.0 or mCVconv != 1.0) {
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        specificHeat[i] *= mCVconv;
        rhoPtr[i] *= mRhoConv;
        Tptr[i] *= mTconv;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LEOS<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {

  // Do the LEOS operation
  LEOSLookup(soundSpeed, massDensity, specificThermalEnergy,
             LEOS_Cs, mMaterialNumber, mDatabaseName,
             mVelConv, mRhoConv, mEconv);
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
LEOS<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {

  const auto n = gamma.size();
  if (n > 0u) {

    // Terrible hack!  We convert the Fields in-place to LEOS units, which
    // means casting away constness.  Boo!
    double *rhoPtr = const_cast<double*>(&massDensity[0]);
    double *epsPtr = const_cast<double*>(&specificThermalEnergy[0]);

    // Convert to LEOS units
    if (mRhoConv != 1.0 or mEconv != 1.0) {
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        rhoPtr[i] /= mRhoConv;
        epsPtr[i] /= mEconv;
      }
    }

    // Get the LEOS functions
    auto fepsPtr = LEOS_bundle::functionPtr(LEOS_Et, mMaterialNumber, mDatabaseName);

    LEOS_LookupOptions lopt;
    lopt.calculateFunction(true).calculateDFDX(true).calculateDFDY(true);

    // Do the LEOS lookup
    std::vector<double> T(n), dedr(n);
    fepsPtr->inverseEval(n, rhoPtr, epsPtr, &T[0], &dedr[0], &gamma[0], &lopt);

    // Convert back to our runtime units and finish calculation of gamma
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      rhoPtr[i] *= mRhoConv;
      epsPtr[i] *= mEconv;
      gamma[i] *= mCVconv;
      const auto nDen = massDensity[i]/mAtomicWeight;
      gamma[i] = 1.0 + mConstants.molarGasConstant()*nDen*safeInvVar(gamma[i]);
    }
  }
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension>
void
LEOS<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {

  const auto n = bulkModulus.size();
  if (n > 0u) {
    Field<Dimension, Scalar> dummy("blago", bulkModulus);  // scratch field
    this->setPressureAndDerivs(dummy, dummy, bulkModulus, massDensity, specificThermalEnergy);
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      bulkModulus[i] = std::abs(bulkModulus[i]*massDensity[i]);
    }
  }
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LEOS<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {

  // Do the LEOS operation
  LEOSLookup(entropy, massDensity, specificThermalEnergy,
             LEOS_St, mMaterialNumber, mDatabaseName,
             mSconv, mRhoConv, mEconv);
}

//------------------------------------------------------------------------------
// Set the melt temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
LEOS<Dimension>::
setMeltTemperature(Field<Dimension, Scalar>& meltTemperature,
                   const Field<Dimension, Scalar>& massDensity,
                   const Field<Dimension, Scalar>& specificThermalEnergy) const {

  // Do the LEOS operation
  LEOSLookup(meltTemperature, massDensity, specificThermalEnergy,
             LEOS_Tm, mMaterialNumber, mDatabaseName,
             mTconv, mRhoConv, mEconv);
}

//------------------------------------------------------------------------------
// Pressure (single value)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LEOS<Dimension>::
pressure(const Scalar massDensity, 
         const Scalar specificThermalEnergy) const {
  return this->applyPressureLimits(LEOSLookupValue(massDensity, specificThermalEnergy,
                                                   LEOS_Pt, mMaterialNumber, mDatabaseName,
                                                   mPconv, mRhoConv, mEconv));
}

//------------------------------------------------------------------------------
// Pressure, DP/Du, DP/Drho (single value)
//------------------------------------------------------------------------------
template<typename Dimension>
std::tuple<typename Dimension::Scalar,
           typename Dimension::Scalar,
           typename Dimension::Scalar>
LEOS<Dimension>::
pressureAndDerivs(const Scalar massDensity, 
                  const Scalar specificThermalEnergy) const {
  double P, DPDu, DPDrho;
  LEOSLookupValue(P, DPDu, DPDrho,
                  massDensity, specificThermalEnergy,
                  LEOS_Pt, mMaterialNumber, mDatabaseName,
                  mPconv, mRhoConv, mEconv);
  const auto Plim = this->applyPressureLimits(P);
  if (P == Plim) {
    return std::make_tuple(P, DPDu, DPDrho);
  } else {
    return std::make_tuple(Plim, 0.0, 0.0);
  }
}

//------------------------------------------------------------------------------
// Temperature (single value)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LEOS<Dimension>::
temperature(const Scalar massDensity, 
            const Scalar specificThermalEnergy) const {

  // Convert to LEOS units
  const L8DOUBLE rho = massDensity/mRhoConv;
  const L8DOUBLE eps = specificThermalEnergy/mEconv;

  // Get the LEOS functions
  auto epsPtr = LEOS_bundle::functionPtr(LEOS_Et, mMaterialNumber, mDatabaseName);

  // Do the LEOS lookup
  auto T = epsPtr->inverseEval(rho, eps);

  // Return in our units
  return mTconv * T;
}

//------------------------------------------------------------------------------
// Specific thermal energy (single value)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LEOS<Dimension>::
specificThermalEnergy(const Scalar massDensity, 
                      const Scalar temperature) const {

  // Convert to LEOS units
  const L8DOUBLE rho = massDensity/mRhoConv;
  const L8DOUBLE T = temperature/mTconv;

  // Get the LEOS functions
  auto epsPtr = LEOS_bundle::functionPtr(LEOS_Et, mMaterialNumber, mDatabaseName);

  // Do the LEOS lookup
  auto eps = epsPtr->eval(rho, T);

  // Return in our units
  return mEconv * eps;
}

//------------------------------------------------------------------------------
// Specific heat (single value)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LEOS<Dimension>::
specificHeat(const Scalar massDensity, 
             const Scalar temperature) const {

  // Convert to LEOS units
  const L8DOUBLE rho = massDensity/mRhoConv;
  const L8DOUBLE T = temperature/mTconv;

  // Get the LEOS functions
  auto epsPtr = LEOS_bundle::functionPtr(LEOS_Et, mMaterialNumber, mDatabaseName);

  LEOS_LookupOptions lopt;
  lopt.calculateFunction(true).calculateDFDX(true).calculateDFDY(true);

  // Do the LEOS lookup
  double eps, DEDT, DEDrho;
  epsPtr->eval(rho, T, eps, DEDT, DEDrho, &lopt);

  // Return in our units
  return DEDT * mCVconv;
}

//------------------------------------------------------------------------------
// Sound speed (single value)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LEOS<Dimension>::
soundSpeed(const Scalar massDensity, 
           const Scalar specificThermalEnergy) const {
  return LEOSLookupValue(massDensity, specificThermalEnergy,
                         LEOS_Cs, mMaterialNumber, mDatabaseName,
                         mVelConv, mRhoConv, mEconv);
}

//------------------------------------------------------------------------------
// Gamma (single value)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LEOS<Dimension>::
gamma(const Scalar massDensity, 
      const Scalar specificThermalEnergy) const {

  // Convert to LEOS units
  const L8DOUBLE rho = massDensity/mRhoConv;
  const L8DOUBLE eps = specificThermalEnergy/mEconv;

  // Get the LEOS functions
  auto epsPtr = LEOS_bundle::functionPtr(LEOS_Et, mMaterialNumber, mDatabaseName);

  LEOS_LookupOptions lopt;
  lopt.calculateFunction(true).calculateDFDX(true).calculateDFDY(true);

  // Do the LEOS lookup
  double T, DEDT, DEDrho;
  epsPtr->inverseEval(rho, eps, T, DEDrho, DEDT, &lopt);

  // Return the result
  DEDT *= mCVconv;
  const auto nDen = massDensity/mAtomicWeight;
  return 1.0 + mConstants.molarGasConstant()*nDen*safeInvVar(DEDT);
}

//------------------------------------------------------------------------------
// Bulk modulus (single value)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LEOS<Dimension>::
bulkModulus(const Scalar massDensity, 
            const Scalar specificThermalEnergy) const {
  auto [P, DPDu, DPDrho] = this->pressureAndDerivs(massDensity, specificThermalEnergy);
  return std::abs(massDensity*DPDrho);  // Protect from negative values
}

//------------------------------------------------------------------------------
// Entropy (single value)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LEOS<Dimension>::
entropy(const Scalar massDensity, 
        const Scalar specificThermalEnergy) const {
  return LEOSLookupValue(massDensity, specificThermalEnergy,
                         LEOS_St, mMaterialNumber, mDatabaseName,
                         mSconv, mRhoConv, mEconv);
}

//------------------------------------------------------------------------------
// Melt temperature (single value)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LEOS<Dimension>::
meltTemperature(const Scalar massDensity, 
                const Scalar specificThermalEnergy) const {

  // For reasons I don't understand looking up the single value of Tm is not
  // working (returns an uninitialized value), but the array interface does.
  vector<double> Tm(1, -1.0), rho(1, massDensity), eps(1, specificThermalEnergy);
  LEOSLookup(Tm, rho, eps,
             LEOS_Tm, mMaterialNumber, mDatabaseName,
             mTconv, mRhoConv, mEconv);
  return Tm[0];

  // return LEOSLookupValue(massDensity, specificThermalEnergy,
  //                        LEOS_Tm, mMaterialNumber, mDatabaseName,
  //                        mTconv, mRhoConv, mEconv);
}

//------------------------------------------------------------------------------
// Do the inverse lookup for eps(rho, P)
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
LEOS<Dimension>::
specificThermalEnergyForPressure(const Scalar Ptarget,
                                 const Scalar rho,
                                 const Scalar epsMin,
                                 const Scalar epsMax,
                                 const Scalar epsTol,
                                 const Scalar Ptol,
                                 const unsigned maxIterations,
                                 const bool verbose) const {

  // Convert to LEOS units
  const L8DOUBLE rho_ = rho/mRhoConv;
  const L8DOUBLE P_   = Ptarget/mPconv;

  // Get the LEOS function
  auto fepsPtr = LEOS_bundle::functionPtr(LEOS_Et, mMaterialNumber, mDatabaseName);
  auto fPptr = LEOS_bundle::functionPtr(LEOS_Pt, mMaterialNumber, mDatabaseName);

  // Do the LEOS lookup
  auto T_ = fPptr->inverseEval(rho_, P_);
  auto eps_ = fepsPtr->eval(rho_, T_);

  // Return in our units
  return mEconv * eps_;
}

//------------------------------------------------------------------------------
// Lookup the reference temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
double
LEOS<Dimension>::referenceTemperature() const {
  L8DOUBLE T0;
  auto matPtr = LEOS_bundle::materialPtr(mMaterialNumber, mDatabaseName);
  matPtr->getMetadata(MAT_TEMP_0, T0);
  return T0 * mTconv;
}

//------------------------------------------------------------------------------
// Return the material descriptor string
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
LEOS<Dimension>::descriptor() const {
  std::ostringstream os;
  auto matPtr = LEOS_bundle::materialPtr(mMaterialNumber, mDatabaseName);
  os << *matPtr;
  return os.str();
}

}
