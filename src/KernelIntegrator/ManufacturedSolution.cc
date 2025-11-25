//---------------------------------Spheral++----------------------------------//
// ManufacturedSolution
//
// Represents some manufactured solutions along with integration coefficients
//----------------------------------------------------------------------------//
#include "ManufacturedSolution.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// ManufacturedFunction
//------------------------------------------------------------------------------
template<typename Dimension>
ManufacturedFunction<Dimension>::
ManufacturedFunction() {
}

template<typename Dimension>
double
ManufacturedFunction<Dimension>::
evaluateCoefficient(const KernelIntegrationData<Dimension>& kid) const {
  return this->evaluate(kid.time, kid.ordinate);
}

//------------------------------------------------------------------------------
// ManufacturedConstantFunction
//------------------------------------------------------------------------------
template<typename Dimension>
ManufacturedConstantFunction<Dimension>::
ManufacturedConstantFunction(const double coefficient):
  mCoefficient(coefficient) {
}

template<typename Dimension>
double
ManufacturedConstantFunction<Dimension>::
evaluate(const double t, const Vector& x) const {
  return mCoefficient;
}

template<typename Dimension>
typename Dimension::Vector
ManufacturedConstantFunction<Dimension>::
evaluateSpatialGradient(const double t, const Vector& x) const {
  return Vector::zero();
}

template<typename Dimension>
typename Dimension::SymTensor
ManufacturedConstantFunction<Dimension>::
evaluateSpatialHessian(const double t, const Vector& x) const {
  return SymTensor::zero();
}

template<typename Dimension>
double
ManufacturedConstantFunction<Dimension>::
evaluateTimeDerivative(const double t, const Vector& x) const {
  return 0.0;
}

//------------------------------------------------------------------------------
// ManufacturedSteadyStateFunction
//------------------------------------------------------------------------------
template<typename Dimension>
ManufacturedSteadyStateFunction<Dimension>::
ManufacturedSteadyStateFunction(const double t,
                                std::shared_ptr<ManufacturedFunction<Dimension>> func):
  mT(t),
  mFunc(func) {
}

template<typename Dimension>
double
ManufacturedSteadyStateFunction<Dimension>::
evaluate(const double t, const Vector& x) const {
  return mFunc->evaluate(mT, x);
}

template<typename Dimension>
typename Dimension::Vector
ManufacturedSteadyStateFunction<Dimension>::
evaluateSpatialGradient(const double t, const Vector& x) const {
  return mFunc->evaluateSpatialGradient(mT, x);
}

template<typename Dimension>
typename Dimension::SymTensor
ManufacturedSteadyStateFunction<Dimension>::
evaluateSpatialHessian(const double t, const Vector& x) const {
  return mFunc->evaluateSpatialHessian(mT, x);
}

template<typename Dimension>
double
ManufacturedSteadyStateFunction<Dimension>::
evaluateTimeDerivative(const double t, const Vector& x) const {
  return 0.0;
}

//------------------------------------------------------------------------------
// ManufacturedSinusoidalFunction
//------------------------------------------------------------------------------
template<typename Dimension>
ManufacturedSinusoidalFunction<Dimension>::
ManufacturedSinusoidalFunction(const std::vector<double>& coefficients):
  mCoefficients(coefficients) {
  VERIFY(mCoefficients.size() == Dimension::nDim + 2);
}

template<typename Dimension>
double
ManufacturedSinusoidalFunction<Dimension>::
evaluate(const double t, const Vector& x) const {
  auto val = 1.0;
  for (auto d = 0; d < Dimension::nDim; ++d) {
    val *= std::cos(mCoefficients[d + 2] * (x[d] + t));
  }
  return mCoefficients[0] * (1 + mCoefficients[1] * val);
}

template<typename Dimension>
typename Dimension::Vector
ManufacturedSinusoidalFunction<Dimension>::
evaluateSpatialGradient(const double t, const Vector& x) const {
  auto val = Vector::one();
  for (auto d1 = 0; d1 < Dimension::nDim; ++d1) {
    for (auto d2 = 0; d2 < Dimension::nDim; ++d2) {
      if (d1 == d2) {
        val[d1] *= -mCoefficients[d1 + 2] * std::sin(mCoefficients[d1 + 2] * (x[d1] + t));
      }
      else {
        val[d1] *= std::cos(mCoefficients[d2 + 2] * (x[d2] + t));
      }
    }
  }
  return mCoefficients[0] * mCoefficients[1] * val;
}

template<typename Dimension>
typename Dimension::SymTensor
ManufacturedSinusoidalFunction<Dimension>::
evaluateSpatialHessian(const double t, const Vector& x) const {
  auto sinVal = Vector::zero();
  auto cosVal = Vector::zero();
  auto cosProd = 1.0;
  for (auto d = 0; d < Dimension::nDim; ++d) {
    const auto arg = mCoefficients[d + 2] * (x[d] + t);
    sinVal[d] = std::sin(arg);
    cosVal[d] = std::cos(arg);
    cosProd *= cosVal[d];
  }

  // Initialize value
  auto val = SymTensor::one();
  
  // Diagonal terms
  for (auto d = 0; d < Dimension::nDim; ++d) {
    val(d, d) = -mCoefficients[d + 2] * mCoefficients[d + 2] * cosProd;
  }

  // Off-diagonal terms
  for (auto d1 = 0; d1 < Dimension::nDim - 1; ++d1) {
    for (auto d2 = d1+1; d2 < Dimension::nDim; ++d2) {
      for (auto d3 = 0; d3 < Dimension::nDim; ++d3) {
        if (d1 == d3 || d2 == d3) {
          val(d1, d2) *= sinVal[d3] * mCoefficients[d3 + 2];;
        }
        else {
          val(d1, d2) *= cosVal[d3];
        }
      }
    }
  }
  
  return mCoefficients[0] * mCoefficients[1] * val;
}

template<typename Dimension>
double
ManufacturedSinusoidalFunction<Dimension>::
evaluateTimeDerivative(const double t, const Vector& x) const {
  const auto valDx = evaluateSpatialGradient(t, x);
  return valDx.sumElements();
}

//------------------------------------------------------------------------------
// ManufacturedWaveFunction
//------------------------------------------------------------------------------
template<typename Dimension>
ManufacturedWaveFunction<Dimension>::
ManufacturedWaveFunction(const std::vector<double>& coefficients):
  mCoefficients(coefficients) {
  VERIFY(mCoefficients.size() == 2);
}

template<typename Dimension>
double
ManufacturedWaveFunction<Dimension>::
evaluate(const double t, const Vector& x) const {
  const auto dist = x.magnitude();
  const auto k0 = mCoefficients[0];
  const auto k1 = mCoefficients[1];
  const auto tdist = dist - t;
  const auto tdist2 = tdist * tdist;
  const auto t2 = t * t;
  const auto t26 = 6. + t2;
  return k0*(1. + 1./(std::exp(k1*tdist2)*t26));
}

template<typename Dimension>
typename Dimension::Vector
ManufacturedWaveFunction<Dimension>::
evaluateSpatialGradient(const double t, const Vector& x) const {
  const auto dist = x.magnitude();
  const auto k0 = mCoefficients[0];
  const auto k1 = mCoefficients[1];
  const auto tdist = dist - t;
  const auto tdist2 = tdist * tdist;
  const auto t2 = t * t;
  const auto t26 = 6. + t2;
  auto val = Vector::zero();
  for (auto d = 0; d < Dimension::nDim; ++d) {
    val[d] = (-2*k0*k1*x[d]*(-t + dist))/(std::exp(k1*tdist2)*t26*dist);
  }
  return val;
}

template<typename Dimension>
typename Dimension::SymTensor
ManufacturedWaveFunction<Dimension>::
evaluateSpatialHessian(const double t, const Vector& x) const {
  const auto dist = x.magnitude();
  const auto dist2 = dist * dist;
  const auto dist3 = dist * dist * dist;
  const auto k0 = mCoefficients[0];
  const auto k1 = mCoefficients[1];
  const auto tdist = dist - t;
  const auto tdist2 = tdist * tdist;
  const auto t2 = t * t;
  const auto t26 = 6. + t2;
  auto val = SymTensor::zero();
  for (auto d1 = 0; d1 < Dimension::nDim; ++d1) {
    for (auto d2 = d1; d2 < Dimension::nDim; ++d2) {
      if (d1 == d2) {
        const auto xd2 = x[d1] * x[d1];
        val(d1, d2) = (2*k0*k1*(t*(dist2 - xd2) - dist3 + 2*k1*xd2*(t2*dist - 2*t*(dist2) + dist3)))/(std::exp(k1*tdist2)*(t26)*dist3);
      }
      else {
        val(d1, d2) = (2*k0*k1*x[d1]*x[d2]*(2*k1*t2*dist + 2*k1*dist3 - t*(1 + 4*k1*(dist2))))/(std::exp(k1*tdist2)*(t26)*dist3);
      }
    }
  }
  return val;
}

template<typename Dimension>
double
ManufacturedWaveFunction<Dimension>::
evaluateTimeDerivative(const double t, const Vector& x) const {
  const auto dist = x.magnitude();
  const auto k0 = mCoefficients[0];
  const auto k1 = mCoefficients[1];
  const auto tdist = dist - t;
  const auto tdist2 = tdist * tdist;
  const auto t2 = t * t;
  const auto t3 = t2 * t;
  const auto t26 = 6. + t2;
  const auto t262 = t26 * t26;
  return (-2*k0*(t + 6*k1*t + k1*t3 - 6*k1*dist - k1*t2*dist))/(std::exp(k1*tdist2)*t262);
}

//------------------------------------------------------------------------------
// ManufacturedTransportSolution
//------------------------------------------------------------------------------
template<typename Dimension>
ManufacturedTransportSolution<Dimension>::
ManufacturedTransportSolution(const double c,
                              const int numOrdinates,
                              const double angularNorm,
                              const std::vector<Vector>& ordinates,
                              std::shared_ptr<ManufacturedFunction<Dimension>> psiFunc,
                              std::shared_ptr<ManufacturedFunction<Dimension>> sigmaAFunc):
  mCInv(1. / c),
  mNumOrdinates(numOrdinates),
  mAngularNorm(angularNorm),
  mOrdinates(ordinates),
  mPsiFunc(psiFunc),
  mSigmaAFunc(sigmaAFunc) {
}

template<typename Dimension>
double
ManufacturedTransportSolution<Dimension>::
evaluatePhi(const double t, const Vector& x) const {
  return mPsiFunc->evaluate(t, x) * mAngularNorm;
}

template<typename Dimension>
std::vector<double>
ManufacturedTransportSolution<Dimension>::
evaluatePsi(const double t, const Vector& x) const {
  const auto psi = mPsiFunc->evaluate(t, x);
  return std::vector<double>(mNumOrdinates, psi);
}

template<typename Dimension>
std::vector<double>
ManufacturedTransportSolution<Dimension>::
evaluateSource(const double t, const Vector& x) const {
  const auto psi = mPsiFunc->evaluate(t, x);
  const auto psiDx = mPsiFunc->evaluateSpatialGradient(t, x);
  const auto psiDt = mPsiFunc->evaluateTimeDerivative(t, x);
  const auto sigmaA = mSigmaAFunc->evaluate(t, x);

  const auto isoSource = mCInv * psiDt + sigmaA * psi;
  std::vector<double> source(mNumOrdinates);
  for (auto o = 0; o < mNumOrdinates; ++o) {
    source[o] = isoSource + mOrdinates[o].dot(psiDx);
  }
  
  return source;
}

template<typename Dimension>
std::vector<double>
ManufacturedTransportSolution<Dimension>::
evaluateCoefficient(const KernelIntegrationData<Dimension>& kid) const {
  return evaluateSource(kid.time, kid.ordinate);
}

} // end namespace Spheral
