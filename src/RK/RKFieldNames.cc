#include "RK/RKFieldNames.hh"

namespace Spheral {

const std::string RKFieldNames::rkOrders = "rkOrders";
const std::string RKFieldNames::rkCorrectionsBase = "rkCorrections/";
const std::string RKFieldNames::reproducingKernelBase = "reproducingKernel/";
const std::string RKFieldNames::rkCorrections(const RKOrder order)     { return RKFieldNames::rkCorrectionsBase + std::to_string(static_cast<int>(order)); }
const std::string RKFieldNames::reproducingKernel(const RKOrder order) { return RKFieldNames::reproducingKernelBase + std::to_string(static_cast<int>(order)); }
RKOrder RKFieldNames::correctionOrder(const std::string& x)      { const auto i = x.find("/"); return static_cast<RKOrder>(std::stoi(x.substr(i+1))); }

}
