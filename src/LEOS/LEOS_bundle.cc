//---------------------------------Spheral++----------------------------------//
// LEOS_bundle
//
// A singleton object which maintains the necessary data for intializing the 
// LEOS package.
//
// Mon May  6 16:18:06 PDT 2013
//----------------------------------------------------------------------------//
#include "LEOS_bundle.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/range.hh"

#include <algorithm>
using std::cout;
using std::cerr;
using std::endl;

namespace Spheral {

//------------------------------------------------------------------------------
// Get the instance.
//------------------------------------------------------------------------------
LEOS_bundle&
LEOS_bundle::
instance() {
  static LEOS_bundle theInstance;
  return theInstance;
}

//------------------------------------------------------------------------------
// Get an LEOS data base
//------------------------------------------------------------------------------
LEOS::LEOS_DatabasePtr
LEOS_bundle::
databasePtr(const name_t dbname,
            const file_t dbformat) {
  LEOS_bundle& bundle = LEOS_bundle::instance();
  if (bundle.mDatabasePtrs.find(dbname) == bundle.mDatabasePtrs.end()) {
    LEOS::LEOS_DatabaseOptions opts;
    opts.publicDatafileType(dbformat);
    auto dbPtr = LEOS::getDatabase(dbname, &opts);
    dbPtr->setup();
    bundle.mDatabasePtrs[dbname] = dbPtr;
    CHECK2(dbPtr->isValid(),
           "Database " << dbname << " invalid : " << dbPtr->getStatus()->message());
  }
  return bundle.mDatabasePtrs[dbname];
}

//------------------------------------------------------------------------------
// Get an LEOS material
//------------------------------------------------------------------------------
LEOS::LEOS_MaterialPtr
LEOS_bundle::
materialPtr(const eosnum_t eosnum,
            const name_t dbname,
            const file_t dbformat) {
  LEOS_bundle& bundle = LEOS_bundle::instance();
  const auto key = std::make_pair(eosnum, dbname);
  if (bundle.mMatPtrs.find(key) == bundle.mMatPtrs.end()) {
    auto dbPtr = bundle.databasePtr(dbname, dbformat);
    auto matPtr = dbPtr->getMaterial(eosnum);
    bundle.mMatPtrs[key] = matPtr;
    
    // Add the EOS functions we need
    LEOS::LEOS_FunctionOptions opts;
    opts.interpolantMemoryType(LEOS::MT_CPU);
    opts.interpolantTemporaryMemoryType(LEOS::MT_CPU);
    opts.execSpace(LEOS::ExecType_t::CPU);
    // opts.setTcalc(true);
    for (const auto& funcName: bundle.mFuncTemplate) {
      auto funcPtr = matPtr->getFunction(funcName, LEOS::BIMOND, &opts);
      if (funcPtr->isValid()) bundle.mFuncPtrs[key][funcName] = funcPtr;
    }
    // CHECK(bundle.mFuncPtrs[key].size() == bundle.mFuncTemplate.size());

    // Initialize the material
    matPtr->initializeFunctions();
    {
      BEGIN_CONTRACT_SCOPE;
      for (auto [funcName, funcPtr]: bundle.mFuncPtrs[key]) {
        CHECK2(funcPtr->isValid(),
               "Function " << funcName << " invalid : " << funcPtr->getStatus()->message());
      }
      END_CONTRACT_SCOPE;
    }
  }
  return bundle.mMatPtrs[key];
}

//------------------------------------------------------------------------------
// Get an LEOS function
//------------------------------------------------------------------------------
LEOS::LEOS_FunctionPtr
LEOS_bundle::
functionPtr(const name_t fname, const eosnum_t eosnum, const name_t dbname) {
  REQUIRE(LEOS::isInitialized());
  LEOS_bundle& bundle = LEOS_bundle::instance();
  const auto key = std::make_pair(eosnum, dbname);
  CHECK(bundle.mFuncPtrs.find(key) != bundle.mFuncPtrs.end());
  auto itr = bundle.mFuncPtrs[key].find(fname);
  VERIFY2(itr != bundle.mFuncPtrs[key].end(),
          "LEOS Error: " << fname << " is not available for EOS " << eosnum << " in " << dbname);
  auto funcPtr = itr->second;
  ENSURE(funcPtr->isValid());
  return funcPtr;
}

//------------------------------------------------------------------------------
// Private constructor.
//------------------------------------------------------------------------------
LEOS_bundle::
LEOS_bundle():
  mFuncTemplate({LEOS::LEOS_Et,      // total specific thermal energy
                 // LEOS::LEOS_Ei,      // ionic energy
                 // LEOS::LEOS_Ee,      // electron energy
                 // LEOS::LEOS_E2p,     // two-phase energy
                 LEOS::LEOS_Pt,      // Total pressure
                 LEOS::LEOS_Tm,      // Melt temperature
                 LEOS::LEOS_Cs,      // Sound speed
                 LEOS::LEOS_St}),    // Entropy
  mDatabasePtrs(),
  mMatPtrs(),
  mFuncPtrs() {

  // Initialize LEOS.
  // We make LEOS work in CGS units internally, which we convert to the user selected
  // units through our interface.
  LEOS::LEOS_StartupOptions opts;
  opts.units(LEOS::LEOS_UNITS_CGS);
#ifdef SPHERAL_ENABLE_MPI
  opts.communicator(Communicator::communicator());
#endif
#ifdef ENABLE_HIP
  opts.permanentMemoryType(LEOS::MT_CPU);
  opts.temporaryMemoryType(LEOS::MT_CPU);
  opts.execSpace(LEOS::ExecType_t::CPU);
#endif
  LEOS::startup(opts);
}

//------------------------------------------------------------------------------
// Private destructor.
//------------------------------------------------------------------------------
LEOS_bundle::
~LEOS_bundle() {
  LEOS::shutdown();  // Tell LEOS to close down and clean up
}

}
