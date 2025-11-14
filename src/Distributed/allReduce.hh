//---------------------------------Spheral++----------------------------------//
// allReduce
//
// Hide (some) of the details about doing MPI all reduces.
//
// Created by JMO, Wed Feb 10 14:38:05 PST 2010
//----------------------------------------------------------------------------//
#ifndef __Spheral_allReduce__
#define __Spheral_allReduce__

#include "Utilities/DataTypeTraits.hh"
#include "Communicator.hh"

#ifdef SPHERAL_ENABLE_MPI
#include <mpi.h>
#endif

namespace Spheral {
#ifdef SPHERAL_ENABLE_MPI
//------------------------------------------------------------------------------
// MPI version
//------------------------------------------------------------------------------

#define SPHERAL_OP_MIN MPI_MIN
#define SPHERAL_OP_MAX MPI_MAX
#define SPHERAL_OP_SUM MPI_SUM
#define SPHERAL_OP_PROD MPI_PROD
#define SPHERAL_OP_LAND MPI_LAND
#define SPHERAL_OP_LOR MPI_LOR
#define SPHERAL_OP_MINLOC MPI_MINLOC
#define SPHERAL_OP_MAXLOC MPI_MAXLOC

template<typename Value>
constexpr Value
allReduce(const Value& value, const MPI_Op op,
          const MPI_Comm comm = Communicator::communicator()) {
  CHECK(!(op == SPHERAL_OP_MINLOC || op == SPHERAL_OP_MAXLOC));
  Value tmp = value;
  Value result;
  MPI_Allreduce(&tmp, &result, 1,
                DataTypeTraits<Value>::MpiDataType(), op, comm);
  return result;
}

template<typename Value>
constexpr std::pair<Value, int>
allReduceLoc(const Value value, const MPI_Op op,
             const MPI_Comm comm = Communicator::communicator()) {
  CHECK(op == SPHERAL_OP_MINLOC || op == SPHERAL_OP_MAXLOC);
  struct {
    Value val;
    int rank;
  } in, out;

  MPI_Comm_rank(comm, &in.rank);
  in.val = value;

  MPI_Allreduce(&in, &out, 1, DataTypeTraits<Value>::MpiLocDataType(), op, comm);

  return {out.val, out.rank};
}


template<typename Value>
constexpr Value
distScan(const Value& value, const MPI_Op op,
     const MPI_Comm comm = Communicator::communicator()) {
  CHECK(!(op == SPHERAL_OP_MINLOC || op == SPHERAL_OP_MAXLOC));
  Value tmp = value;
  Value result;
  MPI_Scan(&tmp, &result, 1, DataTypeTraits<Value>::MpiDataType(), op, comm);
  return result;
}

inline void
Barrier(const MPI_Comm comm = Communicator::communicator()) {
  MPI_Barrier(comm);
}

#else
//------------------------------------------------------------------------------
// Non-MPI version
//------------------------------------------------------------------------------

#define SPHERAL_OP_MIN 1
#define SPHERAL_OP_MAX 2
#define SPHERAL_OP_SUM 3
#define SPHERAL_OP_PROD 4
#define SPHERAL_OP_LAND 5
#define SPHERAL_OP_LOR 6
#define SPHERAL_OP_MINLOC 7
#define SPHERAL_OP_MAXLOC 8

template<typename Value>
constexpr Value
allReduce(const Value& value, const int /*op*/, const int = 0) {
  return value;
}

template<typename Value>
inline std::pair<Value, int>
allReduceLoc(const Value value, const int /*op*/,
             const int = 0) {
  return {value, 0};
}

template<typename Value>
constexpr Value
distScan(const Value& value, const int /*op*/, const int = 0) {
  return value;
}

inline void
Barrier(const int = 0) {
  return;
}
#endif
}
#endif

