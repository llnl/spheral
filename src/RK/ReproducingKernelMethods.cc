#include "RK/ReproducingKernelMethods.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
ReproducingKernelMethods<Dimension>::
ReproducingKernelMethods(const RKOrder order):
  mOrder(order) {
  switch(order) {
  case RKOrder::ZerothOrder:
    mEvaluateBaseKernel = &RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateBaseKernel;
    mEvaluateBaseGradient = &RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateBaseGradient;
    mEvaluateBaseHessian = &RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateBaseHessian;
    mEvaluateBaseKernelAndGradient = &RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateBaseKernelAndGradient;
    mEvaluateKernel = &RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateKernel;
    mEvaluateGradient = &RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateGradient;
    mEvaluateHessian = &RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateHessian;
    mEvaluateKernelAndGradient = &RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateKernelAndGradient;
    mEvaluateKernelAndGradients = &RKUtilities<Dimension, RKOrder::ZerothOrder>::evaluateKernelAndGradients;
    mComputeCorrections = &RKUtilities<Dimension, RKOrder::ZerothOrder>::computeCorrections;
    mComputeNormal = &RKUtilities<Dimension, RKOrder::ZerothOrder>::computeNormal;
    mGetTransformationMatrix = &RKUtilities<Dimension, RKOrder::ZerothOrder>::getTransformationMatrix;
    mApplyTransformation = RKUtilities<Dimension, RKOrder::ZerothOrder>::applyTransformation;
    mGradCorrectionsSize = RKUtilities<Dimension, RKOrder::ZerothOrder>::gradCorrectionsSize;
    mHessCorrectionsSize = RKUtilities<Dimension, RKOrder::ZerothOrder>::hessCorrectionsSize;
    break;
  case RKOrder::LinearOrder:
    mEvaluateBaseKernel = &RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateBaseKernel;
    mEvaluateBaseGradient = &RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateBaseGradient;
    mEvaluateBaseHessian = &RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateBaseHessian;
    mEvaluateBaseKernelAndGradient = &RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateBaseKernelAndGradient;
    mEvaluateKernel = &RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateKernel;
    mEvaluateGradient = &RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateGradient;
    mEvaluateHessian = &RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateHessian;
    mEvaluateKernelAndGradient = &RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateKernelAndGradient;
    mEvaluateKernelAndGradients = &RKUtilities<Dimension, RKOrder::LinearOrder>::evaluateKernelAndGradients;
    mComputeCorrections = &RKUtilities<Dimension, RKOrder::LinearOrder>::computeCorrections;
    mComputeNormal = &RKUtilities<Dimension, RKOrder::LinearOrder>::computeNormal;
    mGetTransformationMatrix = &RKUtilities<Dimension, RKOrder::LinearOrder>::getTransformationMatrix;
    mApplyTransformation = RKUtilities<Dimension, RKOrder::LinearOrder>::applyTransformation;
    mGradCorrectionsSize = RKUtilities<Dimension, RKOrder::LinearOrder>::gradCorrectionsSize;
    mHessCorrectionsSize = RKUtilities<Dimension, RKOrder::LinearOrder>::hessCorrectionsSize;
    break;
  case RKOrder::QuadraticOrder:
    mEvaluateBaseKernel = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateBaseKernel;
    mEvaluateBaseGradient = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateBaseGradient;
    mEvaluateBaseHessian = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateBaseHessian;
    mEvaluateBaseKernelAndGradient = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateBaseKernelAndGradient;
    mEvaluateKernel = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateKernel;
    mEvaluateGradient = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateGradient;
    mEvaluateHessian = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateHessian;
    mEvaluateKernelAndGradient = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateKernelAndGradient;
    mEvaluateKernelAndGradients = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::evaluateKernelAndGradients;
    mComputeCorrections = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::computeCorrections;
    mComputeNormal = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::computeNormal;
    mGetTransformationMatrix = &RKUtilities<Dimension, RKOrder::QuadraticOrder>::getTransformationMatrix;
    mApplyTransformation = RKUtilities<Dimension, RKOrder::QuadraticOrder>::applyTransformation;
    mGradCorrectionsSize = RKUtilities<Dimension, RKOrder::QuadraticOrder>::gradCorrectionsSize;
    mHessCorrectionsSize = RKUtilities<Dimension, RKOrder::QuadraticOrder>::hessCorrectionsSize;
    break;
  case RKOrder::CubicOrder:
    mEvaluateBaseKernel = &RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateBaseKernel;
    mEvaluateBaseGradient = &RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateBaseGradient;
    mEvaluateBaseHessian = &RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateBaseHessian;
    mEvaluateBaseKernelAndGradient = &RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateBaseKernelAndGradient;
    mEvaluateKernel = &RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateKernel;
    mEvaluateGradient = &RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateGradient;
    mEvaluateHessian = &RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateHessian;
    mEvaluateKernelAndGradient = &RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateKernelAndGradient;
    mEvaluateKernelAndGradients = &RKUtilities<Dimension, RKOrder::CubicOrder>::evaluateKernelAndGradients;
    mComputeCorrections = &RKUtilities<Dimension, RKOrder::CubicOrder>::computeCorrections;
    mComputeNormal = &RKUtilities<Dimension, RKOrder::CubicOrder>::computeNormal;
    mGetTransformationMatrix = &RKUtilities<Dimension, RKOrder::CubicOrder>::getTransformationMatrix;
    mApplyTransformation = RKUtilities<Dimension, RKOrder::CubicOrder>::applyTransformation;
    mGradCorrectionsSize = RKUtilities<Dimension, RKOrder::CubicOrder>::gradCorrectionsSize;
    mHessCorrectionsSize = RKUtilities<Dimension, RKOrder::CubicOrder>::hessCorrectionsSize;
    break;
  case RKOrder::QuarticOrder:
    mEvaluateBaseKernel = &RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateBaseKernel;
    mEvaluateBaseGradient = &RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateBaseGradient;
    mEvaluateBaseHessian = &RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateBaseHessian;
    mEvaluateBaseKernelAndGradient = &RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateBaseKernelAndGradient;
    mEvaluateKernel = &RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateKernel;
    mEvaluateGradient = &RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateGradient;
    mEvaluateHessian = &RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateHessian;
    mEvaluateKernelAndGradient = &RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateKernelAndGradient;
    mEvaluateKernelAndGradients = &RKUtilities<Dimension, RKOrder::QuarticOrder>::evaluateKernelAndGradients;
    mComputeCorrections = &RKUtilities<Dimension, RKOrder::QuarticOrder>::computeCorrections;
    mComputeNormal = &RKUtilities<Dimension, RKOrder::QuarticOrder>::computeNormal;
    mGetTransformationMatrix = &RKUtilities<Dimension, RKOrder::QuarticOrder>::getTransformationMatrix;
    mApplyTransformation = RKUtilities<Dimension, RKOrder::QuarticOrder>::applyTransformation;
    mGradCorrectionsSize = RKUtilities<Dimension, RKOrder::QuarticOrder>::gradCorrectionsSize;
    mHessCorrectionsSize = RKUtilities<Dimension, RKOrder::QuarticOrder>::hessCorrectionsSize;
    break;
  case RKOrder::QuinticOrder:
    mEvaluateBaseKernel = &RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateBaseKernel;
    mEvaluateBaseGradient = &RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateBaseGradient;
    mEvaluateBaseHessian = &RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateBaseHessian;
    mEvaluateBaseKernelAndGradient = &RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateBaseKernelAndGradient;
    mEvaluateKernel = &RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateKernel;
    mEvaluateGradient = &RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateGradient;
    mEvaluateHessian = &RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateHessian;
    mEvaluateKernelAndGradient = &RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateKernelAndGradient;
    mEvaluateKernelAndGradients = &RKUtilities<Dimension, RKOrder::QuinticOrder>::evaluateKernelAndGradients;
    mComputeCorrections = &RKUtilities<Dimension, RKOrder::QuinticOrder>::computeCorrections;
    mComputeNormal = &RKUtilities<Dimension, RKOrder::QuinticOrder>::computeNormal;
    mGetTransformationMatrix = &RKUtilities<Dimension, RKOrder::QuinticOrder>::getTransformationMatrix;
    mApplyTransformation = RKUtilities<Dimension, RKOrder::QuinticOrder>::applyTransformation;
    mGradCorrectionsSize = RKUtilities<Dimension, RKOrder::QuinticOrder>::gradCorrectionsSize;
    mHessCorrectionsSize = RKUtilities<Dimension, RKOrder::QuinticOrder>::hessCorrectionsSize;
    break;
  case RKOrder::SexticOrder:
    mEvaluateBaseKernel = &RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateBaseKernel;
    mEvaluateBaseGradient = &RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateBaseGradient;
    mEvaluateBaseHessian = &RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateBaseHessian;
    mEvaluateBaseKernelAndGradient = &RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateBaseKernelAndGradient;
    mEvaluateKernel = &RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateKernel;
    mEvaluateGradient = &RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateGradient;
    mEvaluateHessian = &RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateHessian;
    mEvaluateKernelAndGradient = &RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateKernelAndGradient;
    mEvaluateKernelAndGradients = &RKUtilities<Dimension, RKOrder::SexticOrder>::evaluateKernelAndGradients;
    mComputeCorrections = &RKUtilities<Dimension, RKOrder::SexticOrder>::computeCorrections;
    mComputeNormal = &RKUtilities<Dimension, RKOrder::SexticOrder>::computeNormal;
    mGetTransformationMatrix = &RKUtilities<Dimension, RKOrder::SexticOrder>::getTransformationMatrix;
    mApplyTransformation = RKUtilities<Dimension, RKOrder::SexticOrder>::applyTransformation;
    mGradCorrectionsSize = RKUtilities<Dimension, RKOrder::SexticOrder>::gradCorrectionsSize;
    mHessCorrectionsSize = RKUtilities<Dimension, RKOrder::SexticOrder>::hessCorrectionsSize;
    break;
  case RKOrder::SepticOrder:
    mEvaluateBaseKernel = &RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateBaseKernel;
    mEvaluateBaseGradient = &RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateBaseGradient;
    mEvaluateBaseHessian = &RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateBaseHessian;
    mEvaluateBaseKernelAndGradient = &RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateBaseKernelAndGradient;
    mEvaluateKernel = &RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateKernel;
    mEvaluateGradient = &RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateGradient;
    mEvaluateHessian = &RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateHessian;
    mEvaluateKernelAndGradient = &RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateKernelAndGradient;
    mEvaluateKernelAndGradients = &RKUtilities<Dimension, RKOrder::SepticOrder>::evaluateKernelAndGradients;
    mComputeCorrections = &RKUtilities<Dimension, RKOrder::SepticOrder>::computeCorrections;
    mComputeNormal = &RKUtilities<Dimension, RKOrder::SepticOrder>::computeNormal;
    mGetTransformationMatrix = &RKUtilities<Dimension, RKOrder::SepticOrder>::getTransformationMatrix;
    mApplyTransformation = RKUtilities<Dimension, RKOrder::SepticOrder>::applyTransformation;
    mGradCorrectionsSize = RKUtilities<Dimension, RKOrder::SepticOrder>::gradCorrectionsSize;
    mHessCorrectionsSize = RKUtilities<Dimension, RKOrder::SepticOrder>::hessCorrectionsSize;
    break;
  default:
    VERIFY2("Unknown order passed to ReproducingKernelMethods", false);
  }
}

//------------------------------------------------------------------------------
// Default constructor (everything NULL, invalid)
//------------------------------------------------------------------------------
template<typename Dimension>
ReproducingKernelMethods<Dimension>::
ReproducingKernelMethods():
  mOrder(RKOrder::ZerothOrder),
  mGradCorrectionsSize(-1),
  mHessCorrectionsSize(-1),
  mEvaluateBaseKernel(nullptr),
  mEvaluateBaseGradient(nullptr),
  mEvaluateBaseHessian(nullptr),
  mEvaluateBaseKernelAndGradient(nullptr),
  mEvaluateKernel(nullptr),
  mEvaluateGradient(nullptr),
  mEvaluateHessian(nullptr),
  mEvaluateKernelAndGradient(nullptr),
  mEvaluateKernelAndGradients(nullptr),
  mComputeCorrections(nullptr),
  mComputeNormal(nullptr),
  mGetTransformationMatrix(nullptr),
  mApplyTransformation(nullptr) {
}

//------------------------------------------------------------------------------
// Copy
//------------------------------------------------------------------------------
template<typename Dimension>
ReproducingKernelMethods<Dimension>::
ReproducingKernelMethods(const ReproducingKernelMethods<Dimension>& rhs):
  mOrder(rhs.mOrder),
  mGradCorrectionsSize(rhs.mGradCorrectionsSize),
  mHessCorrectionsSize(rhs.mHessCorrectionsSize),
  mEvaluateBaseKernel(rhs.mEvaluateBaseKernel),
  mEvaluateBaseGradient(rhs.mEvaluateBaseGradient),
  mEvaluateBaseHessian(rhs.mEvaluateBaseHessian),
  mEvaluateBaseKernelAndGradient(rhs.mEvaluateBaseKernelAndGradient),
  mEvaluateKernel(rhs.mEvaluateKernel),
  mEvaluateGradient(rhs.mEvaluateGradient),
  mEvaluateHessian(rhs.mEvaluateHessian),
  mEvaluateKernelAndGradient(rhs.mEvaluateKernelAndGradient),
  mEvaluateKernelAndGradients(rhs.mEvaluateKernelAndGradients),
  mComputeCorrections(rhs.mComputeCorrections),
  mComputeNormal(rhs.mComputeNormal),
  mGetTransformationMatrix(rhs.mGetTransformationMatrix),
  mApplyTransformation(rhs.mApplyTransformation) {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
ReproducingKernelMethods<Dimension>&
ReproducingKernelMethods<Dimension>::
operator=(const ReproducingKernelMethods<Dimension>& rhs) {
  mOrder = rhs.mOrder;
  mGradCorrectionsSize = rhs.mGradCorrectionsSize;
  mHessCorrectionsSize = rhs.mHessCorrectionsSize;
  mEvaluateBaseKernel = rhs.mEvaluateBaseKernel;
  mEvaluateBaseGradient = rhs.mEvaluateBaseGradient;
  mEvaluateBaseHessian = rhs.mEvaluateBaseHessian;
  mEvaluateBaseKernelAndGradient = rhs.mEvaluateBaseKernelAndGradient;
  mEvaluateKernel = rhs.mEvaluateKernel;
  mEvaluateGradient = rhs.mEvaluateGradient;
  mEvaluateHessian = rhs.mEvaluateHessian;
  mEvaluateKernelAndGradient = rhs.mEvaluateKernelAndGradient;
  mEvaluateKernelAndGradients = rhs.mEvaluateKernelAndGradients;
  mComputeCorrections = rhs.mComputeCorrections;
  mComputeNormal = rhs.mComputeNormal;
  mGetTransformationMatrix = rhs.mGetTransformationMatrix;
  mApplyTransformation = rhs.mApplyTransformation;
  return *this;
}

//------------------------------------------------------------------------------
// Equivalence
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ReproducingKernelMethods<Dimension>::
operator==(const ReproducingKernelMethods<Dimension>& rhs) const {
  return (mOrder == rhs.mOrder and
          mGradCorrectionsSize == rhs.mGradCorrectionsSize and
          mHessCorrectionsSize == rhs.mHessCorrectionsSize);
}

}
