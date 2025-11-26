#include "DataTypeFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// ConstantIntegrationCoefficient
//------------------------------------------------------------------------------
template<typename Dimension, typename CoefficientType>
inline
ConstantIntegrationCoefficient<Dimension, CoefficientType>::
ConstantIntegrationCoefficient():
  mDataSet(false) {
}

template<typename Dimension, typename CoefficientType>
inline
ConstantIntegrationCoefficient<Dimension, CoefficientType>::
ConstantIntegrationCoefficient(CoefficientType coeff):
  mDataSet(true),
  mCoeff(coeff) {
}

template<typename Dimension, typename CoefficientType>
inline
CoefficientType
ConstantIntegrationCoefficient<Dimension, CoefficientType>::
evaluateCoefficient(const KernelIntegrationData<Dimension>& /* kid */) const {
  return mCoeff;
}

template<typename Dimension, typename CoefficientType>
inline
void
ConstantIntegrationCoefficient<Dimension, CoefficientType>::
setData(CoefficientType coeff) {
  mDataSet = true;
  mCoeff = coeff;
}

template<typename Dimension, typename CoefficientType>
inline
const CoefficientType&
ConstantIntegrationCoefficient<Dimension, CoefficientType>::
getData() const {
  CHECK(mDataSet);
  return mCoeff;
}

//------------------------------------------------------------------------------
// DefaultIntegrationCoefficient
//------------------------------------------------------------------------------
template<typename Dimension, typename CoefficientType>
inline
CoefficientType
DefaultIntegrationCoefficient<Dimension, CoefficientType>::
evaluateCoefficient(const KernelIntegrationData<Dimension>& /* kid */) const {
  return CoefficientType::one();
}

template<>
inline
Dim<1>::Scalar
DefaultIntegrationCoefficient<Dim<1>, Dim<1>::Scalar>::
evaluateCoefficient(const KernelIntegrationData<Dim<1>>& /* kid */) const {
  return 1.0;
}
template<>
inline
Dim<2>::Scalar
DefaultIntegrationCoefficient<Dim<2>, Dim<2>::Scalar>::
evaluateCoefficient(const KernelIntegrationData<Dim<2>>& /* kid */) const {
  return 1.0;
}
template<>
inline
Dim<3>::Scalar
DefaultIntegrationCoefficient<Dim<3>, Dim<3>::Scalar>::
evaluateCoefficient(const KernelIntegrationData<Dim<3>>& /* kid */) const {
  return 1.0;
}

template<>
inline
std::vector<double>
DefaultIntegrationCoefficient<Dim<1>, std::vector<double>>::
evaluateCoefficient(const KernelIntegrationData<Dim<1>>& /* kid */) const {
  return std::vector<double>();
}
template<>
inline
std::vector<double>
DefaultIntegrationCoefficient<Dim<2>, std::vector<double>>::
evaluateCoefficient(const KernelIntegrationData<Dim<2>>& /* kid */) const {
  return std::vector<double>();
}
template<>
inline
std::vector<double>
DefaultIntegrationCoefficient<Dim<3>, std::vector<double>>::
evaluateCoefficient(const KernelIntegrationData<Dim<3>>& /* kid */) const {
  return std::vector<double>();
}

template<>
inline
std::vector<Dim<1>::Vector>
DefaultIntegrationCoefficient<Dim<1>, std::vector<Dim<1>::Vector>>::
evaluateCoefficient(const KernelIntegrationData<Dim<1>>& /* kid */) const {
  return std::vector<Dim<1>::Vector>();
}
template<>
inline
std::vector<Dim<2>::Vector>
DefaultIntegrationCoefficient<Dim<2>, std::vector<Dim<2>::Vector>>::
evaluateCoefficient(const KernelIntegrationData<Dim<2>>& /* kid */) const {
  return std::vector<Dim<2>::Vector>();
}
template<>
inline
std::vector<Dim<3>::Vector>
DefaultIntegrationCoefficient<Dim<3>, std::vector<Dim<3>::Vector>>::
evaluateCoefficient(const KernelIntegrationData<Dim<3>>& /* kid */) const {
  return std::vector<Dim<3>::Vector>();
}

//------------------------------------------------------------------------------
// FieldListIntegrationCoefficient
// f\left(x\right)=\sum_{i}V_{i}f_{i}W_{i}\left(x\right)
//------------------------------------------------------------------------------
template<typename Dimension, typename CoefficientType>
inline
FieldListIntegrationCoefficient<Dimension, CoefficientType>::
FieldListIntegrationCoefficient() :
  IntegrationCoefficient<Dimension, CoefficientType>(),
  mDataSet(false),
  mData() {
}

template<typename Dimension, typename CoefficientType>
inline
FieldListIntegrationCoefficient<Dimension, CoefficientType>::
FieldListIntegrationCoefficient(const FieldList<Dimension, CoefficientType>& data):
  IntegrationCoefficient<Dimension, CoefficientType>(),
  mDataSet (true),
  mData(data) {
}

template<typename Dimension, typename CoefficientType>
inline
const FieldList<Dimension, CoefficientType>& 
FieldListIntegrationCoefficient<Dimension, CoefficientType>::
getData() const {
  CHECK(mDataSet);
  return mData;
}

template<typename Dimension, typename CoefficientType>
inline
void
FieldListIntegrationCoefficient<Dimension, CoefficientType>::
setData(const FieldList<Dimension, CoefficientType>& data) {
  mData = data;
  mDataSet = true;
}

template<typename Dimension, typename CoefficientType>
inline
void
FieldListIntegrationCoefficient<Dimension, CoefficientType>::
setDataPoint(int nodeListi, int nodei, const CoefficientType& data) {
  CHECK(mDataSet);
  mData(nodeListi, nodei) = data;
}

template<typename Dimension, typename CoefficientType>
inline
CoefficientType
FieldListIntegrationCoefficient<Dimension, CoefficientType>::
evaluateCoefficient(const KernelIntegrationData<Dimension>& kid) const {
  CHECK(mDataSet);
  const auto numIndices = kid.indices.size();
  CHECK(kid.volume.size() == numIndices);
  CHECK(kid.values.size() == numIndices);
  auto val = mData(kid.nodeIndex0.first, kid.nodeIndex0.second); // example value for us to zero out
  DataTypeFunctions<CoefficientType>::zeroOutData(val);
  for (auto i = 0u; i < numIndices; ++i) {
    const auto nodeListi = kid.nodeIndices[i].first;
    const auto nodei = kid.nodeIndices[i].second;
    auto localVal = mData(nodeListi, nodei);
    CHECK(DataTypeFunctions<CoefficientType>::size(val) == DataTypeFunctions<CoefficientType>::size(localVal));
    DataTypeFunctions<CoefficientType>::multiplyData(localVal, kid.volume[i] * kid.values[i]);
    DataTypeFunctions<CoefficientType>::addToData(val, localVal);
  }
  return val;
}

} // end namespace Spheral
