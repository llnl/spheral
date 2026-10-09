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

#include <algorithm>
#include <type_traits>
#include <vector>

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

namespace AllReduceDetail {
// Is Value a non-arithmetic type built from doubles (Vector, Tensor, ...)?
template<typename Value, typename = void>
struct hasDoubleElements: std::false_type {};

template<typename Value>
struct hasDoubleElements<Value, std::void_t<typename DataTypeTraits<Value>::ElementType>>:
    std::bool_constant<not std::is_arithmetic<Value>::value and
                       std::is_same<typename DataTypeTraits<Value>::ElementType, double>::value> {};
}

template<typename Value>
Value
allReduce(const Value& value, const MPI_Op op,
          const MPI_Comm comm = Communicator::communicator()) {
  CHECK(!(op == SPHERAL_OP_MINLOC || op == SPHERAL_OP_MAXLOC));

  // The predefined MPI reduction operations are only defined for the basic
  // MPI datatypes, not the derived types we register for our geometric types
  // (Vector, Tensor, etc.).  MPI rejects those reductions (silently, under
  // MPI_ERRORS_RETURN), so we handle them here: sums are done elementwise,
  // and min/max use the type's own comparison, consistent with the local
  // Field::localMin/localMax.
  if constexpr (AllReduceDetail::hasDoubleElements<Value>::value) {
    const int n = DataTypeTraits<Value>::numElements(value);
    CHECK(sizeof(Value) == n*sizeof(double));
    if (op == SPHERAL_OP_SUM) {
      Value result(value);
      MPI_Allreduce(MPI_IN_PLACE, &result, n, MPI_DOUBLE, op, comm);
      return result;
    } else if (op == SPHERAL_OP_MIN or op == SPHERAL_OP_MAX) {
      int nprocs;
      MPI_Comm_size(comm, &nprocs);
      std::vector<Value> values(nprocs);
      Value tmp = value;
      MPI_Allgather(&tmp, n, MPI_DOUBLE, &values.front(), n, MPI_DOUBLE, comm);
      return (op == SPHERAL_OP_MIN ?
              *std::min_element(values.begin(), values.end()) :
              *std::max_element(values.begin(), values.end()));
    }
  }

  Value tmp = value;
  Value result;
  MPI_Allreduce(&tmp, &result, 1,
                DataTypeTraits<Value>::MpiDataType(), op, comm);
  return result;
}

template<typename Value>
std::pair<Value, int>
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
Value
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
Value
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
Value
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
