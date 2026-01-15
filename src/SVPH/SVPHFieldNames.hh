//---------------------------------Spheral++----------------------------------//
// SVPHFieldNames -- A collection of standard Field names for the SVPH 
// physics package.
//
// Created by JMO, Wed Jan  8 11:15:42 PST 2020
//----------------------------------------------------------------------------//
#ifndef _Spheral_SVPHFieldNames_
#define _Spheral_SVPHFieldNames_

#include <string>

namespace Spheral {

struct SVPHFieldNames {
  const inline static std::string A_SVPH = "A SVPH correction";
  const inline static std::string B_SVPH = "B SVPH correction";
  const inline static std::string gradB_SVPH = "gradB SVPH correction";
};

}

#endif
