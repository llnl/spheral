namespace Spheral {

//------------------------------------------------------------------------------
// Return empty statistics
//------------------------------------------------------------------------------
inline
std::vector<std::shared_ptr<IncrementalStatistic<double>>>
LinearSolver::
statistics() const {
  return std::vector<std::shared_ptr<IncrementalStatistic<double>>>();
}

//------------------------------------------------------------------------------
// Set/get description
//------------------------------------------------------------------------------
inline
void
LinearSolver::
setDescription(std::string desc) {
  mDescription = desc;
}

inline
std::string
LinearSolver::
getDescription() const {
  return mDescription;
}

} // end namespace Spheral
