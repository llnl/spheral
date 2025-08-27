namespace Spheral {

//------------------------------------------------------------------------------
// Return statistics for number of iterations and tolerance reached
//------------------------------------------------------------------------------
inline
std::vector<std::shared_ptr<IncrementalStatistic<double>>>
HypreLinearSolver::
statistics() const {
  return {mIterationStatistics, mFinalResidualStatistics};
}

} // end namespace Spheral

