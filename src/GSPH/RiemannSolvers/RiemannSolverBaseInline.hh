namespace Spheral {

//------------------------------------------------------------------------------
// Get slope limiter
//------------------------------------------------------------------------------
template<typename Dimension>
inline
LimiterBase<Dimension>&
RiemannSolverBase<Dimension>::
limiter() const {
  return mSlopeLimiter;
}

//------------------------------------------------------------------------------
// Get wave speed
//------------------------------------------------------------------------------
template<typename Dimension>
inline
WaveSpeedBase<Dimension>&
RiemannSolverBase<Dimension>::
waveSpeed() const {
  return mWaveSpeed;
}

//------------------------------------------------------------------------------
// set/get linear reconstruction switch
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
RiemannSolverBase<Dimension>::
linearReconstruction() const {
  return mLinearReconstruction;
}

template<typename Dimension>
inline
void
RiemannSolverBase<Dimension>::
linearReconstruction(bool x) {
  mLinearReconstruction=x;
}

}