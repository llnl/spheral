//---------------------------------Spheral++----------------------------------//
// IncrementalStatistic
//
// Allows tracking of incrementally added stats, such as iteration counts
// For now, all variables use the same DataType
// The mean guess prevents overflow for data with a large mean
//----------------------------------------------------------------------------//
#ifndef __Spheral_IncrementalStatistic_hh__
#define __Spheral_IncrementalStatistic_hh__

#include <string>

namespace Spheral {

template<typename DataType>
class IncrementalStatistic {
public:
  // Constructor
  IncrementalStatistic(const DataType meanGuess = 0,
                       const std::string name = "IncrementalStatistic",
                       bool printAdd = false);

  // Add a point of data
  void add(const DataType data);

  // Calculated data
  DataType mean() const;
  DataType variance() const;
  DataType total() const;
  
  // Stored data
  std::string name() const;
  DataType meanGuess() const;
  DataType min() const;
  DataType max() const;
  DataType numPoints() const;
  DataType shiftedSum() const;
  DataType shiftedSum2() const;
  
  // Print to console
  void print() const;

  // Print summary to console
  void printSummary() const;

  // Set name
  void setName(std::string newName);

private:
  // Convenience function
  void getStats(DataType& min,
                DataType& max,
                DataType& numPoints,
                DataType& mean,
                DataType& variance,
                DataType& total) const;
  
  // Input data
  const DataType mMeanGuess;
  std::string mName;
  bool mPrint;
  
  // Incrementally updated data
  DataType mMin;
  DataType mMax;
  DataType mNumPoints;
  DataType mShiftedSum;
  DataType mShiftedSum2;
};

} // end namespace Spheral

#include "IncrementalStatisticInline.hh"

#endif
