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
bmLim(const Scalar rho,
      const Scalar bm) const {
  // bm = rho cs^2, so clamp to ensure cs in [csMin, csMax]
  const auto bmMin = std::max(0.0, rho * mCsMin * mCsMin);
  const auto bmMax = std::min(std::numeric_limits<double>::max(), rho * mCsMax * mCsMax);
  return std::min(std::max(bm, bmMin), bmMax);
}

template<typename Dimension>
inline
double
SingularityEquationOfState<Dimension>::
cvLim(const Scalar cv) const {
  // With this ordering, NaN should return mCvMin
  return std::min(std::max(cv, mCvMin), mCvMax);
}


}
