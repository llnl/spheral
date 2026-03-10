//---------------------------------Spheral++----------------------------------//
// RestartRegistrar
// A singleton object which maintains a list of the restartable objects in
// Spheral.
//
// Created by JMO, Wed May 27 13:44:32 PDT 2009
//----------------------------------------------------------------------------//
#include "RestartRegistrar.hh"
#include "FileIO/FileIO.hh"
#include "Utilities/removeElements.hh"

#include <algorithm>
#include <functional>
#include <iostream>
#include <sstream>

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Get the instance.
//------------------------------------------------------------------------------
RestartRegistrar&
RestartRegistrar::
instance() {
  static RestartRegistrar theInstance;
  return theInstance;
}

//------------------------------------------------------------------------------
// Register a RestartHandle.
//------------------------------------------------------------------------------
void
RestartRegistrar::
registerRestartHandle(std::shared_ptr<RestartHandle> restartHandlePtr,
                      const unsigned priority) {
  this->removeExpiredPointers();
  CHECK(mPriorities.size() == mRestartHandles.size());
  std::weak_ptr<RestartHandle> wptr(restartHandlePtr);
  if (not haveRestartHandle(restartHandlePtr)) {
    priority_iterator itr = upper_bound(mPriorities.begin(), mPriorities.end(), priority);
    const size_t delta = distance(mPriorities.begin(), itr);
    mRestartHandles.insert(mRestartHandles.begin() + delta, wptr);
    mPriorities.insert(itr, priority);
  }
  ENSURE(haveRestartHandle(restartHandlePtr));
  ENSURE(mRestartHandles.size() == mPriorities.size());
}

//------------------------------------------------------------------------------
// Unregister a RestartHandle.
//------------------------------------------------------------------------------
void
RestartRegistrar::
unregisterRestartHandle(std::shared_ptr<RestartHandle> restartHandlePtr) {
  this->removeExpiredPointers();
  std::weak_ptr<RestartHandle> wptr(restartHandlePtr);
  VERIFY(haveRestartHandle(restartHandlePtr));
  iterator itr = find_if(this->begin(), this->end(),
                         [wptr] (const std::weak_ptr<RestartHandle> val) {
                           // Weak pointers don't have operator==, so we have to provide something.
                           return val.lock() == wptr.lock();
                         });
  CHECK(itr != this->end());
  const size_t delta = distance(this->begin(), itr);
  mRestartHandles.erase(itr);
  mPriorities.erase(mPriorities.begin() + delta);
  ENSURE(not haveRestartHandle(restartHandlePtr));
  ENSURE(mRestartHandles.size() == mPriorities.size());
}

//------------------------------------------------------------------------------
// Check whether the given RestartHandle is registered.
//------------------------------------------------------------------------------
bool
RestartRegistrar::
haveRestartHandle(const std::shared_ptr<RestartHandle> restartHandlePtr) const {
  // const_iterator itr = std::find_if(this->begin(), this->end(),
  //                                   [restartHandlePtr] (const std::weak_ptr<RestartHandle> val) {
  //                                     // Weak pointers don't have operator==, so we have to provide something.
  //                                     return val.lock() == restartHandlePtr.lock();
  //                                   });
  // return (itr != this->end());
  const_iterator itr = this->begin();
  while (itr < this->end() and
         itr->lock() != restartHandlePtr) ++itr;
  return (itr != this->end());
}

//------------------------------------------------------------------------------
// Eliminate any pointers to expired objects.
//------------------------------------------------------------------------------
void
RestartRegistrar::
removeExpiredPointers() {
  vector<size_t> expiredIndicies;
  for (size_t i = 0; i != mRestartHandles.size(); ++i) {
    if (mRestartHandles[i].expired()) expiredIndicies.push_back(i);
  }
  removeElements(mRestartHandles, expiredIndicies);
  removeElements(mPriorities, expiredIndicies);

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    ENSURE(mRestartHandles.size() == mPriorities.size());
    for (const_iterator itr = this->begin();
         itr != this->end();
         ++itr ) ENSURE(not itr->expired());
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Generate unique labels for each of the objects we currently are restarting.
//------------------------------------------------------------------------------
vector<string>
RestartRegistrar::
uniqueLabels() {
  this->removeExpiredPointers();
  vector<string> result;
  unsigned counter = 0;
  for (const_iterator itr = this->begin();
       itr != this->end();
       ++itr) {
    string l = itr->lock()->label();
    if (find(result.begin(), result.end(), l) != result.end()) {
      std::stringstream newlabel;
      newlabel << l << "_" << counter;
      l = newlabel.str();
      ++counter;
    }
    CHECK(find(result.begin(), result.end(), l) == result.end());
    result.push_back(l);
  }
  ENSURE(result.size() == mRestartHandles.size());
  return result;
}

//------------------------------------------------------------------------------
// Print out the labels.
//------------------------------------------------------------------------------
void
RestartRegistrar::
printLabels() {
  this->removeExpiredPointers();
  const vector<string> labels = this->uniqueLabels();
  for (std::vector<std::string>::const_iterator itr = labels.begin();
       itr != labels.end();
       ++itr) {
    std::cout << *itr << std::endl;
  }
}

//------------------------------------------------------------------------------
// Dump the state of all objects.
//------------------------------------------------------------------------------
void
RestartRegistrar::
dumpState(FileIO& file) {
  this->removeExpiredPointers();
  const vector<string> labels = this->uniqueLabels();
  CHECK(labels.size() == mRestartHandles.size());
  for (size_t i = 0; i != labels.size(); ++i) {
    mRestartHandles[i].lock()->dumpState(file, labels[i]);
  }
}

//------------------------------------------------------------------------------
// Restore the state of all objects.
//------------------------------------------------------------------------------
void
RestartRegistrar::
restoreState(const FileIO& file) {
  this->removeExpiredPointers();
  const vector<string> labels = this->uniqueLabels();
  CHECK(labels.size() == mRestartHandles.size());
  for (auto i = 0u; i != labels.size(); ++i) {

    // We do a bit of chicanery here for backwards-comparability with the Pybindgen version of Spheral.
    auto label = labels[i];
    auto dirPath = file.groupName(label);
    auto varName = file.variableName(label);
    if (varName == "SolidSPHHydroBase" and file.pathExists(dirPath + "/SolidSPHHydroBase_1")) {
      label = dirPath + "/SolidSPHHydroBase_1";
    } else if (varName == "SolidSPHHydroBaseRZ" and file.pathExists(dirPath + "/SolidSPHHydroBaseRZ_1")) {
      label = dirPath + "/SolidSPHHydroBaseRZ_1";
    } else if (varName == "SolidCRKSPHHydroBase" and file.pathExists(dirPath + "SolidCRKSPHHydroBase_1")) {
      label = dirPath + "/SolidCRKSPHHydroBase_1";
    } else if (varName == "SolidCRKSPHHydroBaseRZ" and file.pathExists(dirPath + "SolidCRKSPHHydroBaseRZ_1")) {
      label = dirPath + "/SolidCRKSPHHydroBaseRZ_1";
    }

    mRestartHandles[i].lock()->restoreState(file, label);
  }
}

//------------------------------------------------------------------------------
// Default constructor (private).
//------------------------------------------------------------------------------
RestartRegistrar::
RestartRegistrar():
   mRestartHandles(),
   mPriorities() {
}

//------------------------------------------------------------------------------
// Destructor (private).
//------------------------------------------------------------------------------
RestartRegistrar::
~RestartRegistrar() {
}

}
