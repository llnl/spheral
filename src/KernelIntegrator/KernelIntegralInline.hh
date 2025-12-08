#include "Utilities/DataTypeTraits.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// LinearIntegral
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
LinearIntegral<Dimension, DataType>::
initialize(const FlatConnectivity<Dimension>& flatConnectivity) {
  mValues.assign(flatConnectivity.numNodes(), DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// LinearKernelStdVector
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
LinearKernelStdVector<Dimension>::
initialize(const FlatConnectivity<Dimension>& flatConnectivity) {
  std::vector<double> zero(mSize, 0.0);
  this->mValues.assign(flatConnectivity.numNodes(), zero);
}

//------------------------------------------------------------------------------
// LinearKernelStdVector
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
LinearGradStdVector<Dimension>::
initialize(const FlatConnectivity<Dimension>& flatConnectivity) {
  std::vector<typename Dimension::Vector> zero(mSize, Dimension::Vector::zero());
  this->mValues.assign(flatConnectivity.numNodes(), zero);
}

//------------------------------------------------------------------------------
// BilinearIntegral
//------------------------------------------------------------------------------
template<typename Dimension, typename BaseDataType>
inline
void
BilinearIntegral<Dimension, BaseDataType>::
initialize(const FlatConnectivity<Dimension>& flatConnectivity) {
  const auto numNodes = flatConnectivity.numNodes();
  mValues.resize(numNodes);
  const auto zero = DataTypeTraits<BaseDataType>::zero();
  if (this->volume()) {
    // If there are volume integrals (regardless of surface integrals), we need all points
    for (auto i = 0; i < numNodes; ++i) {
#if REPLACEOVERLAP
      const auto numValues = flatConnectivity.numNeighbors(i);
#else
      const auto numValues = flatConnectivity.numOverlapNeighbors(i);
#endif
      mValues[i].assign(numValues, zero);
    }
  }
  else if (this->surface()) {
    // If there are no volume integrals but there are surface integrals, we only need the surface storage
    for (auto i = 0; i < numNodes; ++i) {
#if REPLACEOVERLAP
      const auto numValues = flatConnectivity.numNeighbors(i);
#else
      const auto numValues = flatConnectivity.numOverlapNeighbors(i);
#endif
      const auto numSurfaces = flatConnectivity.numSurfaces(i);
      if (numSurfaces > 0) {
        mValues[i].assign(numValues, zero);
      }
      else {
        mValues[i].resize(0);
      }
    }
  }
  else {
    // We need something to do
    VERIFY2(false, "need either surface or volume integral terms");
  }
}

//------------------------------------------------------------------------------
// LinearSurfaceDependentIntegral
//------------------------------------------------------------------------------
template<typename Dimension, typename BaseDataType>
inline
void
LinearSurfaceDependentIntegral<Dimension, BaseDataType>::
initialize(const FlatConnectivity<Dimension>& flatConnectivity) {
  const auto numNodes = flatConnectivity.numNodes();
  mValues.resize(numNodes);
  const auto zero = DataTypeTraits<BaseDataType>::zero();
  for (auto i = 0; i < numNodes; ++i) {
    const auto numSurfaces = flatConnectivity.numSurfaces(i);
    if (numSurfaces > 0) {
      mValues[i].assign(numSurfaces, zero);
    }
    else {
      mValues[i].resize(0);
    }
  }
}

//------------------------------------------------------------------------------
// LinearSurfaceNormalKernelStdVector
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
LinearSurfaceNormalKernelStdVector<Dimension>::
initialize(const FlatConnectivity<Dimension>& flatConnectivity) {
  std::vector<typename Dimension::Vector> zero(mSize, Dimension::Vector::zero());
  const auto numNodes = flatConnectivity.numNodes();
  this->mValues.resize(numNodes);
  for (auto i = 0; i < numNodes; ++i) {
    const auto numSurfaces = flatConnectivity.numSurfaces(i);
    if (numSurfaces > 0) {
      this->mValues[i].assign(numSurfaces, zero);
    }
    else {
      this->mValues[i].resize(0);
    }
  }
}

//------------------------------------------------------------------------------
// BilinearSurfaceDependentIntegral
//------------------------------------------------------------------------------
template<typename Dimension, typename BaseDataType>
inline
void
BilinearSurfaceDependentIntegral<Dimension, BaseDataType>::
initialize(const FlatConnectivity<Dimension>& flatConnectivity) {
  const auto numNodes = flatConnectivity.numNodes();
  mValues.resize(numNodes);
  const auto zero = DataTypeTraits<BaseDataType>::zero();
  for (auto i = 0; i < numNodes; ++i) {
#if REPLACEOVERLAP
      const auto numValues = flatConnectivity.numNeighbors(i);
#else
      const auto numValues = flatConnectivity.numOverlapNeighbors(i);
#endif
    const auto numSurfaces = flatConnectivity.numSurfaces(i);
    if (numSurfaces > 0) {
      mValues[i].assign(numValues * numSurfaces, zero);
    }
    else {
      mValues[i].resize(0);
    }
  }
}

//------------------------------------------------------------------------------
// BilinearMultiplyByFieldList
//------------------------------------------------------------------------------
template<typename Dimension, typename BaseDataType>
inline
void
BilinearMultiplyByFieldList<Dimension, BaseDataType>::
finalize(const FlatConnectivity<Dimension>& flatConnectivity) {
  // Copy the base values
  const auto baseVals = mBaseIntegral->values();
  this->mValues = baseVals;

  // Get the FieldList
  const auto mult = this->getCoefficient()->getData();
  
  // For each value, multiply it by the FieldList
  const auto numNodes = flatConnectivity.numNodes();
  for (auto i = 0; i < numNodes; ++i) {
    const auto pairi = flatConnectivity.localToNode(i);
    const auto nodeListi = pairi.first;
    const auto nodei = pairi.second;
    const auto multi = mult(nodeListi, nodei);
    auto& valsi = this->mValues[i];
    const auto numVals = valsi.size();
    for (auto j = 0u; j < numVals; ++j) {
      valsi[j] *= multi;
    }
  }
}

} // end namespace Spheral
