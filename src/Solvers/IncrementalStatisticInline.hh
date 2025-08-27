#include <iostream>
#include <iomanip>
#include <limits>
#include <vector>

#include "Distributed/Communicator.hh"
#include "Utilities/DataTypeTraits.hh"

namespace { // anonymous

template<typename DataType>
inline
DataType
calculateMean(DataType numPoints,
              DataType shiftedSum,
              DataType meanGuess) {
  return shiftedSum / numPoints + meanGuess;
}

template<typename DataType>
inline
DataType
calculateVariance(DataType numPoints,
                  DataType shiftedSum,
                  DataType shiftedSum2) {
  return (shiftedSum2 - shiftedSum * shiftedSum / numPoints) / numPoints;
}

template<typename DataType>
inline
DataType
calculateTotal(DataType numPoints,
               DataType shiftedSum,
               DataType meanGuess) {
  return shiftedSum + numPoints * meanGuess;
}

} // namespace anonymous

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename DataType>
inline
IncrementalStatistic<DataType>::
IncrementalStatistic(const DataType meanGuess,
                     std::string name,
                     bool print):
  mMeanGuess(meanGuess),
  mName(name),
  mPrint(print),
  mMin(std::numeric_limits<DataType>::max()),
  mMax(std::numeric_limits<DataType>::min()),
  mNumPoints(0),
  mShiftedSum(0),
  mShiftedSum2(0) {
}

//------------------------------------------------------------------------------
// Add a data point
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
IncrementalStatistic<DataType>::
add(const DataType data) {
  if (data < mMin) mMin = data;
  if (data > mMax) mMax = data;
  mNumPoints += 1;
  const DataType delta = data - mMeanGuess;
  mShiftedSum += delta;
  mShiftedSum2 += delta * delta;
  
  if (mPrint) {
    std::cout << mName;
    std::cout << "\tadd: " << data;
    std::cout << std::endl;
  }
}

//------------------------------------------------------------------------------
// Calculate mean
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
mean() const {
  return calculateMean(mNumPoints, mShiftedSum, mMeanGuess);
}

//------------------------------------------------------------------------------
// Calculate variance
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
variance() const {
  return calculateVariance(mNumPoints, mShiftedSum, mShiftedSum2);
}

//------------------------------------------------------------------------------
// Calculate total
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
total() const {
  return calculateTotal(mNumPoints, mShiftedSum, mMeanGuess);
}


//------------------------------------------------------------------------------
// Set name
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
IncrementalStatistic<DataType>::
setName(std::string name) {
  mName = name;
  return;
}

//------------------------------------------------------------------------------
// Get mean guess
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
min() const {
  return mMin;
}

//------------------------------------------------------------------------------
// Get mean guess
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
max() const {
  return mMax;
}

//------------------------------------------------------------------------------
// Get mean guess
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
meanGuess() const {
  return mMeanGuess;
}

//------------------------------------------------------------------------------
// Get numPoints
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
numPoints() const {
  return mNumPoints;
}

//------------------------------------------------------------------------------
// Get shiftedSum
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
shiftedSum() const {
  return mShiftedSum;
}

//------------------------------------------------------------------------------
// Get shiftedSum2
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType
IncrementalStatistic<DataType>::
shiftedSum2() const {
  return mShiftedSum2;
}

//------------------------------------------------------------------------------
// Print results
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
IncrementalStatistic<DataType>::
print() const {
  typedef std::pair<std::string, DataType> PairType;

  auto outputPair = [](PairType val) {
    std::cout << std::setw(20) << val.first;
    std::cout << std::setw(14) << std::setprecision(6) << val.second;
  };

  DataType min, max, numPoints, mean, variance, total;
  getStats(min, max, numPoints, mean, variance, total);

  if (Process::getRank() == 0) {
    std::cout << mName << std::endl;
    if (numPoints == 0) {
      PairType val = {"points", numPoints};
      outputPair(val);
      std::cout << "(no data recorded)" << std::endl;
    }
    else {
      // Output in two columns
      std::vector<PairType> vals =
        {{"points",   numPoints}, {"total",   total},
         {"mean",     mean},      {"minimum", min},
         {"variance", variance},  {"maximum", max}};
      auto index = 0;
      for (auto val : vals) {
        outputPair(val);
        if (index % 2 == 1) {
          std::cout << std::endl;
        }
        index += 1;
      }
    }
  }
}

template<typename DataType>
inline
void
IncrementalStatistic<DataType>::
printSummary() const {
  DataType min, max, numPoints, mean, variance, total;
  getStats(min, max, numPoints, mean, variance, total);

  if (Process::getRank() == 0) {
    std::cout << std::left << std::setw(22) << mName;
    if (numPoints == 0) {
      std::cout << " (no data recorded)" << std::endl;
    }
    else {
      std::cout << std::scientific << std::setprecision(4);
      std::cout << " mean=" << mean
                << "  min=" << min
                << "  max=" << max
                << "  tot=" << total << std::endl;
    }
  }
}

// Get raw stats communicated over all ranks
template<typename DataType>
inline
void
IncrementalStatistic<DataType>::
getStats(DataType& min,
         DataType& max,
         DataType& numPoints,
         DataType& mean,
         DataType& variance,
         DataType& total) const {
  // Communicate raw data
  DataType shiftedSum, shiftedSum2;
  const auto mpiDataType = DataTypeTraits<DataType>::MpiDataType();
  MPI_Allreduce(&mMin, &min, 1, mpiDataType, MPI_MIN, Communicator::communicator());
  MPI_Allreduce(&mMax, &max, 1, mpiDataType, MPI_MAX, Communicator::communicator());
  MPI_Allreduce(&mNumPoints, &numPoints, 1, mpiDataType, MPI_SUM, Communicator::communicator());
  MPI_Allreduce(&mShiftedSum, &shiftedSum, 1, mpiDataType, MPI_SUM, Communicator::communicator());
  MPI_Allreduce(&mShiftedSum2, &shiftedSum2, 1, mpiDataType, MPI_SUM, Communicator::communicator());

  // Calculate mean and variance with total num points, since they are divided through
  mean = calculateMean(numPoints, shiftedSum, mMeanGuess);
  variance = calculateVariance(numPoints, shiftedSum, shiftedSum2);
  total = calculateTotal(numPoints, shiftedSum, mMeanGuess);

  // Divide total and numPoints by procs, since these are usually added similarly on every proc
  const auto numProcs = Process::getTotalNumberOfProcesses();
  numPoints /= numProcs;
  total /= numProcs;
}
                        
  


} // end namespace Spheral
