//---------------------------------Spheral++----------------------------------//
// LEOS_bundle
//
// A singleton object which maintains the necessary data for intializing the 
// LEOS package.
//
// Because LEOS works natively in (rho,T) space, but Spheral's EOS interface defaults
// to (rho, eps), we need to prepare for LEOS inverse lookups to get T(rho,eps).
// For this reason we lay out the data in the following orders:
//
// Mon May  6 16:18:06 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_LEOS_bundle__
#define __Spheral_LEOS_bundle__

#include "Utilities/DBC.hh"

#include <LEOS.h>

#include <boost/container_hash/hash.hpp>  // So we can use unordered_map with std::pair as a key

#include <vector>
#include <utility>
#include <unordered_map>

namespace Spheral {

class LEOS_bundle {

public:
  //------------------------===== Public Interface =====-----------------------//
  using name_t = ::LEOS::L8STRING;
  using eosnum_t = ::LEOS::materialType_t;
  using file_t = ::LEOS::fileType;

  // Get the instance.
  static LEOS_bundle& instance();

  // Get an LEOS Database
  static ::LEOS::LEOS_DatabasePtr databasePtr(const name_t dbname = "leos",
                                              const file_t dbformat = ::LEOS::fileType::PDBFILE_LEOS);

  // Get an LEOS material
  static ::LEOS::LEOS_MaterialPtr materialPtr(const eosnum_t eosnum,
                                              const name_t dbname = "leos",
                                              const file_t dbformat = ::LEOS::fileType::PDBFILE_LEOS);

  // Get an LEOS function
  static ::LEOS::LEOS_FunctionPtr functionPtr(const name_t fname,
                                              const eosnum_t eosnum,
                                              const name_t dbname = "leos");

  // No copying or assignment
  LEOS_bundle(const LEOS_bundle&) = delete;
  LEOS_bundle& operator=(const LEOS_bundle&) = delete;

private:
  //------------------------===== Private Interface =====----------------------//
  std::vector<::LEOS::L8STRING>                                                               mFuncTemplate;
  std::unordered_map<name_t, ::LEOS::LEOS_DatabasePtr>                                        mDatabasePtrs;
  std::unordered_map<std::pair<eosnum_t, name_t>, ::LEOS::LEOS_MaterialPtr,
                     boost::hash<std::pair<eosnum_t, name_t>>>                                mMatPtrs;
  std::unordered_map<std::pair<eosnum_t, name_t>, std::map<name_t, ::LEOS::LEOS_FunctionPtr>,
                     boost::hash<std::pair<eosnum_t, name_t>>>                                mFuncPtrs;

  // Private constructor and destructor
  LEOS_bundle();
  ~LEOS_bundle();
};

}

#endif
