//---------------------------------Spheral++----------------------------------//
// RedistributionRegistrar
// A singleton object which maintains a list of methods to notify when a 
// redistribution of the nodes has occurred.
//
// Created by JMO, Thu Dec 10 14:07:01 PST 2009
//----------------------------------------------------------------------------//
#include "RedistributionRegistrar.hh"

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
// Register a RedistributionNotificationHandle.
//------------------------------------------------------------------------------
void
RedistributionRegistrar::
registerRedistributionNotificationHandle(std::shared_ptr<RedistributionNotificationHandle> redistributionHandlePtr) {
  if (not haveRedistributionNotificationHandle(std::weak_ptr<RedistributionNotificationHandle>(redistributionHandlePtr))) {
    mRedistributionNotificationHandles.push_back(std::weak_ptr<RedistributionNotificationHandle>(redistributionHandlePtr));
  }
  ENSURE(haveRedistributionNotificationHandle(std::weak_ptr<RedistributionNotificationHandle>(redistributionHandlePtr)));
}

//------------------------------------------------------------------------------
// Unregister a RedistributionNotificationHandle.
//------------------------------------------------------------------------------
void
RedistributionRegistrar::
unregisterRedistributionNotificationHandle(std::shared_ptr<RedistributionNotificationHandle> redistributionHandlePtr) {
  std::weak_ptr<RedistributionNotificationHandle> wptr(redistributionHandlePtr);
  VERIFY(haveRedistributionNotificationHandle(wptr));
  iterator itr = std::find_if(this->begin(), this->end(),
                              [wptr] (const std::weak_ptr<RedistributionNotificationHandle> val) {
                                // Weak pointers don't have operator==, so we have to provide something.
                                return val.lock() == wptr.lock();
                              });
  CHECK(itr != this->end());
  mRedistributionNotificationHandles.erase(itr);
  ENSURE(not haveRedistributionNotificationHandle(wptr));
}

//------------------------------------------------------------------------------
// Check whether the given RedistributionNotificationHandle is registered.
//------------------------------------------------------------------------------
bool
RedistributionRegistrar::
haveRedistributionNotificationHandle(std::weak_ptr<RedistributionNotificationHandle> redistributionHandlePtr) const {
  const_iterator itr = std::find_if(this->begin(), this->end(),
                                    [redistributionHandlePtr] (const std::weak_ptr<RedistributionNotificationHandle> val) {
                                      // Weak pointers don't have operator==, so we have to provide something.
                                      return val.lock() == redistributionHandlePtr.lock();
                                    });
  return itr != mRedistributionNotificationHandles.end();
}

//------------------------------------------------------------------------------
// Eliminate any pointers to expired objects.
//------------------------------------------------------------------------------
void
RedistributionRegistrar::
removeExpiredPointers() {
  iterator itr = end();
  int i = mRedistributionNotificationHandles.size();
  while (i != 0) {
    --i;
    if (mRedistributionNotificationHandles[i].expired()) {
      mRedistributionNotificationHandles.erase(this->begin() + i);
    }
  }

  CONTRACT_VAR(itr);
  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  {
    for (const_iterator itr = this->begin();
         itr != this->end();
         ++itr ) ENSURE(not itr->expired());
  }
  END_CONTRACT_SCOPE
}

//------------------------------------------------------------------------------
// Notify everyone registered that a redistribution is about to occur.
//------------------------------------------------------------------------------
void
RedistributionRegistrar::
preRedistributionNotifications() {
  removeExpiredPointers();
  for (const_iterator itr = begin(); itr != end(); ++itr) itr->lock()->notifyBeforeRedistribution();
}

//------------------------------------------------------------------------------
// Notify everyone registered that a redistribution has occurred.
//------------------------------------------------------------------------------
void
RedistributionRegistrar::
broadcastRedistributionNotifications() {
  removeExpiredPointers();
  for (const_iterator itr = begin(); itr != end(); ++itr) itr->lock()->notifyAfterRedistribution();
}

//------------------------------------------------------------------------------
// Default constructor (private).
//------------------------------------------------------------------------------
RedistributionRegistrar::
RedistributionRegistrar():
  mRedistributionNotificationHandles() {
}

//------------------------------------------------------------------------------
// Destructor (private).
//------------------------------------------------------------------------------
RedistributionRegistrar::
~RedistributionRegistrar() {
}

}
