//---------------------------------Spheral++----------------------------------//
// PairwiseFieldViewBase -- abstract handle to obfuscate the template types
// for PairwiseFieldView
//
// Created by J. Michael Owen, Wed Jul  8 13:25:21 PDT 2026
//----------------------------------------------------------------------------//
#ifndef _Spheral_NeighborSpace_PairwiseFieldViewBase_
#define _Spheral_NeighborSpace_PairwiseFieldViewBase_

#include "Utilities/GPUUtils.hh"

namespace Spheral {

class PairwiseFieldViewBase:
    public chai::CHAICopyable {

public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors, destructors
  SPHERAL_HOST_DEVICE PairwiseFieldViewBase()                                            = default;
  SPHERAL_HOST_DEVICE PairwiseFieldViewBase(const PairwiseFieldViewBase& rhs)            = default;
  SPHERAL_HOST_DEVICE PairwiseFieldViewBase(PairwiseFieldViewBase&& rhs)                 = default;
  SPHERAL_HOST_DEVICE virtual ~PairwiseFieldViewBase()                                   = default;

  SPHERAL_HOST virtual void move(chai::ExecutionSpace space) = 0;
  SPHERAL_HOST virtual void touch(chai::ExecutionSpace space) = 0;
};

}

#endif

