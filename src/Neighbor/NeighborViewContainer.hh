//---------------------------------Spheral++----------------------------------//
// Common container interface used by NeighborView implementations.
//
// Created by J. Michael Owen, Mon Sep 28 16:53:54 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_NeighborViewContainer__
#define __Spheral_NeighborViewContainer__

#include "chai/ManagedArray.hpp"
#include "chai/ExecutionSpaces.hpp"

#ifdef SPHERAL_ENABLE_GPU
#include <span>
#else
#include <vector>
#endif

namespace Spheral {

#ifdef SPHERAL_ENABLE_GPU

// GPU friendly container for views  
#ifdef SPHERAL_UNIFIED_MEMORY
template<typename T> using NeighborViewContainer = std::span<T>;
#else
template<typename T> using NeighborViewContainer = chai::ManagedArray<T>;
#endif

#else
// Non-GPU container. NeighborView methods can grow this container as needed.
template<typename T> using NeighborViewContainer = std::vector<T>;

#endif

}

#endif
