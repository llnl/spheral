namespace Spheral {

template<typename Dimension>
inline
double
SingularityEquationOfState<Dimension>::
rhoLim(const Scalar rho) const {
  return std::max(mRhoMin, std::min(rho, mRhoMax));
}

template<typename Dimension>
inline
double
SingularityEquationOfState<Dimension>::
eLim(const Scalar rho,
     const Scalar e,
     const SingularityLimitLevel level) const {
  if (mLevel < level) {
    return e;
  }
  auto it = std::lower_bound(mDensitiesCGS.begin(), mDensitiesCGS.end(), rho);
  const auto i = it == mDensitiesCGS.end() ? 0 : std::distance(mDensitiesCGS.begin(), it);
  return std::max(mMinEnergiesCGS[i], e);
}

template<typename Dimension>
inline
double
SingularityEquationOfState<Dimension>::
tLim(const Scalar rho,
     const Scalar t,
     const SingularityLimitLevel level) const {
  if (mLevel < level) {
    return t;
  }
  auto it = std::lower_bound(mDensitiesCGS.begin(), mDensitiesCGS.end(), rho);
  const auto i = it == mDensitiesCGS.end() ? 0 : std::distance(mDensitiesCGS.begin(), it);
  return std::max(mMinTemperaturesCGS[i], t);
}

template<typename Dimension>
inline
double
SingularityEquationOfState<Dimension>::
csLim(const Scalar cs,
      const SingularityLimitLevel level) const {
  if (mLevel < level) {
    return cs;
  }
  // With this ordering, NaN should return mCsMin
  return std::min(std::max(cs, mCsMin), mCsMax);
}

template<typename Dimension>
inline
double
SingularityEquationOfState<Dimension>::
bmLim(const Scalar rho,
      const Scalar bm,
      const SingularityLimitLevel level) const {
  if (mLevel < level) {
    return bm;
  }
  // With this ordering, NaN should return mCsMin
  // bm = rho cs^2
  const auto bmMin = std::max(0.0, rho * mCsMin * mCsMin);
  const auto bmMax = std::min(std::numeric_limits<double>::max(), rho * mCsMax * mCsMax);
  return std::min(std::max(bm, bmMin), bmMax);
}


}
