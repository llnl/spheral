//---------------------------------Spheral++----------------------------------//
// hashes
//
// A collection of simple useful hashes.
//
// Created by JMO, Fri Jan  7 21:56:22 PST 2011
//----------------------------------------------------------------------------//
#ifndef __Spheral_hashes__
#define __Spheral_hashes__

namespace std {

template<typename first, typename second>
struct hash<std::pair<first, second>> {
  size_t operator()(const std::pair<first, second>& p) const {
    return hash<first>()(p.first) ^ hash<second>()(p.second);
  }
};

} // namespace std

#endif

