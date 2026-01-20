//---------------------------------Spheral++----------------------------------//
// RKFieldNames -- A collection of standard Field names for the Reproducing
// Kernel (RK) package.
//
// Created by JMO, Sat Dec 28 13:18:36 PST 2019
//----------------------------------------------------------------------------//
#ifndef _Spheral_RKFieldNames_
#define _Spheral_RKFieldNames_

#include "RK/RKCorrectionParams.hh"
#include <string>

namespace Spheral {

struct RKFieldNames {
  const inline static std::string rkOrders = "rkOrders";
  const inline static std::string rkCorrectionsBase = "rkCorrections_";
  const inline static std::string reproducingKernelBase = "reproducingKernel_";
  const inline static std::string rkCorrections(const RKOrder order)     { return RKFieldNames::rkCorrectionsBase + std::to_string(static_cast<int>(order)); }
  const inline static std::string reproducingKernel(const RKOrder order) { return RKFieldNames::reproducingKernelBase + std::to_string(static_cast<int>(order)); }
  inline static RKOrder correctionOrder(const std::string& x)      { const auto i = x.find("_"); return static_cast<RKOrder>(std::stoi(x.substr(i+1))); }
};

}

#endif
