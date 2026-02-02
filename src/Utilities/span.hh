//---------------------------------Spheral++----------------------------------//
// Wrapper for std::span
//
// If we're using a C++20 compiler we can use std::span, otherwise use
// boost::span
//----------------------------------------------------------------------------//
#ifndef __Spheral_Span__
#define __Spheral_Span__

// Check for C++20 support
// #if __cplusplus >= 202002L
#ifdef SPHERAL_USE_STD_SPAN
  #define SPHERAL_SPAN_HEADER <span>
  #define SPHERAL_SPAN_TYPE std::span
#else
  #define SPHERAL_SPAN_HEADER <boost/core/span.hpp>
  #define SPHERAL_SPAN_TYPE boost::span
#endif

// Include the appropriate header
#include SPHERAL_SPAN_HEADER

#endif
