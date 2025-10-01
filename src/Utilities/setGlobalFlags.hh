//---------------------------------Spheral++----------------------------------//
// setGlobalFlags
//
// Catch all function invoked by SpheralController at problem startup.
//
// Created by JMO, Thu Oct 27 13:23:24 PDT 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_setGlobalFlags__
#define __Spheral_setGlobalFlags__

#ifdef __GNUC__
#include <fenv.h>
#endif

namespace Spheral {
void setGlobalFlags() {

#ifdef __GNUC__
#ifdef SPHERAL_ENABLE_NAN_EXCEPTIONS
  feenableexcept(FE_DIVBYZERO | FE_INVALID);
  // feenableexcept(FE_DIVBYZERO | FE_UNDERFLOW | FE_OVERFLOW | FE_INVALID);
#endif
#endif

}
}
  
#endif
