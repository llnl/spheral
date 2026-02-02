//---------------------------------Spheral++----------------------------------//
// BilinearIndex
//
// Indexes a bilinear form
//----------------------------------------------------------------------------//
#ifndef __Spheral_BilinearIndex_hh__
#define __Spheral_BilinearIndex_hh__

#include <vector>
#include <unordered_map>
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"

namespace Spheral {

// These hashes are specialized for low numbers of integers
// The array version assumes the last index is either 0 or 1
template<typename DataType> struct BilinearHash {
  int operator()(const DataType& v) const;
};

template<typename Dimension>
class BilinearIndex {
public:
  
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Typedefs for bilinear indexing
  typedef typename std::pair<int, int> Pair;
  typedef typename std::unordered_map<Pair, int, BilinearHash<Pair>> MapPair;
  typedef typename std::vector<Pair> VectorPair;
  typedef FieldList<Dimension, MapPair> PairToFlat;
  typedef FieldList<Dimension, VectorPair> FlatToPair;

  // Typedefs for surface indexing
  typedef typename std::array<int, Dimension::nDim> ArrayDim;
  typedef typename std::unordered_map<ArrayDim, int, BilinearHash<ArrayDim>> MapNormal;
  typedef typename std::vector<Vector> NormalType;
  typedef FieldList<Dimension, MapNormal> NormalToFlat;
  typedef FieldList<Dimension, NormalType> FlatToNormal;
  typedef FieldList<Dimension, std::vector<int>> SurfaceFlags;
  
  // Construtor
  BilinearIndex();
  
  // Compute the connectivity for the bilinear form
  void computeConnectivity(const DataBase<Dimension>& dataBase);

  // Return whether connectivity has been computed
  bool connectivityComputed() { return mConnectivityComputed; }

  // Compute the surface indexing
  void computeSurfaceIndexing(const DataBase<Dimension>& dataBase,
                              const State<Dimension>& state);

  bool surfaceIndexingComputed() { return mSurfaceIndexingComputed; }

  // Return the number of nonzero elements on a given row
  int numElementsForRow(const std::pair<int, int>& ni) const;
  
  // For the row ni, for the index nj, get the flat index j
  int flatIndex(const std::pair<int, int>& ni,
                const std::pair<int, int>& nj) const;

  // For the row ni, for the flat index j, get the index nj
  std::pair<int, int> nodeIndex(const std::pair<int, int>& ni,
                                const int& j) const;

  // Get all the indices for this row, helpful for evaluating functions
  const std::vector<std::pair<int, int>>& nodeIndices(const std::pair<int, int>& ni) const;

  // For the row ni, how many surfaces are registered?
  int numSurfacesForRow(const std::pair<int, int>& ni) const;
  
  // For the row ni, get the surface index s with the given normal
  int surfaceIndex(const std::pair<int, int>& ni,
                   const Vector& normal) const;

  // Same as above, but using the key for better efficiency
  int surfaceIndex(const std::pair<int, int>& ni,
                   const ArrayDim& values) const;

  // For the row ni, get the normal directions for each surface 
  const std::vector<Vector>& normal(const std::pair<int, int>& ni) const;

  // Get the facets for this cell that we have flagged for integration
  const std::vector<int>& surfaceFlags(const std::pair<int, int>& ni) const;

  // Convert a normal to an array of ints
  void normalToArray(const Vector& normal,
                     ArrayDim& values) const;
  void arrayToNormal(const ArrayDim& values,
                     Vector& normal) const;
  
private:
  
  // Connectivity for bilinear form
  bool mConnectivityInitialized;
  bool mConnectivityComputed;
  bool mSurfaceIndexingInitialized;
  bool mSurfaceIndexingComputed;
  PairToFlat mBilinearFlatIndex;
  FlatToPair mBilinearNodeIndex;
  FlatToNormal mSurfaceNormal;
  NormalToFlat mSurfaceFlatIndex;
  SurfaceFlags mSurfaceFlags;

  // All normals within this tolerance will be combined for indexing and integrals
  static constexpr double mRoundValue = 1.e8;
  static constexpr double mRoundValueInv = 1.e-8;
  
  // Scratch
  mutable ArrayDim mScratchArray;
}; // end class BilinearIndex

} // end namespace Spheral

#include "BilinearIndexInline.hh"

#endif
