//---------------------------------Spheral++----------------------------------//
// sampleFieldList2LatticeCIC
//
// Use CIC (cloud in cell) interpolation to sample a FieldList to a lattice.
// The results are returned as a 
// vector<vector<Value> >.
//
// Created by JMO, Wed Feb  3 11:02:20 PST 2010
//----------------------------------------------------------------------------//
#include "Geometry/Dimension.hh"
#include "binFieldList2Lattice.hh"
#include "Kernel/TableKernel.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "NodeList/NodeList.hh"
#include "Geometry/MathTraits.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Distributed/allReduce.hh"
#include "Distributed/BoundingVolumeDistributedBoundary.hh"
#include "Distributed/Communicator.hh"

#include "Utilities/DBC.hh"

#include <algorithm>
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

#ifdef SPHERAL_ENABLE_MPI
extern "C" {
#include <mpi.h>
}
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Helper method specialized per dimension to increment grid cell values.
//------------------------------------------------------------------------------
template<typename Value>
void
incrementCellValues(vector<Value>& values, 
                    const Value& value,
                    const Dim<1>::Vector& xi,
                    const Dim<1>::SymTensor& Hi,
                    const Dim<1>::Vector& xmin,
                    const Dim<1>::Vector& xmax,
                    const vector<unsigned>& nsample,
                    const TableKernel<Dim<1> >& W) {
  REQUIRE(nsample.size() == 1 and nsample[0] > 0);
  typedef Dim<1>::Vector Vector;
  const Vector extent = Neighbor<Dim<1> >::HExtent(Hi, W.kernelExtent());
  const double xstep = (xmax.x() - xmin.x())/nsample[0];
  const size_t nx = size_t(extent.x()/xstep);
  const int ipx = max(0, min(int(nsample[0]) - 1, int((xi.x() - xmin.x())/xstep)));
  const double Hdeti = Hi.Determinant();
  for (int idx = -nx; idx != (int)nx + 1; ++idx) {
    const int ix = ipx + idx;
    if (ix > 0 and ix < (int)nsample[0]) {
      const double etaMag = (Hi*Vector(idx*xstep)).magnitude();
      const unsigned icell = ix;
      CHECK(icell < nsample[0]);
      values[icell] += value*W(etaMag, Hdeti)*Hdeti;
    }
  }
}

template<typename Value>
void
incrementCellValues(vector<Value>& values, 
                    const Value& value,
                    const Dim<2>::Vector& xi,
                    const Dim<2>::SymTensor& Hi,
                    const Dim<2>::Vector& xmin,
                    const Dim<2>::Vector& xmax,
                    const vector<unsigned>& nsample,
                    const TableKernel<Dim<2> >& W) {
  REQUIRE(nsample.size() == 2 and nsample[0] > 0 and nsample[1] > 0);
  typedef Dim<2>::Vector Vector;
  const unsigned ncells = nsample[0]*nsample[1];
  const Vector extent = Neighbor<Dim<2> >::HExtent(Hi, W.kernelExtent());
  const double xstep = (xmax.x() - xmin.x())/nsample[0];
  const double ystep = (xmax.y() - xmin.y())/nsample[1];
  const size_t nx = size_t(extent.x()/xstep);
  const size_t ny = size_t(extent.y()/ystep);
  const int ipx = max(0, min(int(nsample[0]) - 1, int((xi.x() - xmin.x())/xstep)));
  const int ipy = max(0, min(int(nsample[1]) - 1, int((xi.y() - xmin.y())/ystep)));
  const double Hdeti = Hi.Determinant();
  for (int idy = -ny; idy != (int)ny + 1; ++idy) {
    const int iy = ipy + idy;
    if (iy >= 0 and iy < (int)nsample[1]) {
      const unsigned iyoff = iy*nsample[0];
      for (int idx = -nx; idx != (int)nx + 1; ++idx) {
        const int ix = ipx + idx;
        if (ix >= 0 and ix < (int)nsample[0]) {
          const double etaMag = (Hi*Vector(idx*xstep, idy*ystep)).magnitude();
          const unsigned icell = ix + iyoff;
          CONTRACT_VAR(ncells);
          CHECK(icell < ncells);
          values[icell] += value*W(etaMag, Hdeti)*Hdeti;
        }
      }
    }
  }
}

template<typename Value>
void
incrementCellValues(vector<Value>& values, 
                    const Value& value,
                    const Dim<3>::Vector& xi,
                    const Dim<3>::SymTensor& Hi,
                    const Dim<3>::Vector& xmin,
                    const Dim<3>::Vector& xmax,
                    const vector<unsigned>& nsample,
                    const TableKernel<Dim<3> >& W) {
  REQUIRE(nsample.size() == 2 and nsample[0] > 0 and nsample[1] > 0 and nsample[2] > 0);
  typedef Dim<3>::Vector Vector;
  const unsigned ncells = nsample[0]*nsample[1]*nsample[2];
  const Vector extent = Neighbor<Dim<3> >::HExtent(Hi, W.kernelExtent());
  const double xstep = (xmax.x() - xmin.x())/nsample[0];
  const double ystep = (xmax.y() - xmin.y())/nsample[1];
  const double zstep = (xmax.z() - xmin.z())/nsample[2];
  const size_t nx = size_t(extent.x()/xstep);
  const size_t ny = size_t(extent.y()/ystep);
  const size_t nz = size_t(extent.z()/zstep);
  const int ipx = max(0, min(int(nsample[0]) - 1, int((xi.x() - xmin.x())/xstep)));
  const int ipy = max(0, min(int(nsample[1]) - 1, int((xi.y() - xmin.y())/ystep)));
  const int ipz = max(0, min(int(nsample[2]) - 1, int((xi.z() - xmin.z())/zstep)));
  const double Hdeti = Hi.Determinant();
  for (int idz = -nz; idz != (int)nz + 1; ++idz) {
    const int iz = ipz + idz;
    if (iz >= 0 and iz < (int)nsample[2]) {
      const unsigned izoff = iz*nsample[0]*nsample[1];
      for (int idy = -ny; idy != (int)ny + 1; ++idy) {
        const int iy = ipy + idy;
        if (iy >= 0 and iy < (int)nsample[1]) {
          const unsigned iyoff = iy*nsample[0];
          for (int idx = -nx; idx != (int)nx + 1; ++idx) {
            const int ix = ipx + idx;
            if (ix > 0 and ix < (int)nsample[0]) {
              const double etaMag = (Hi*Vector(idx*xstep, idy*ystep)).magnitude();
              const unsigned icell = ix + iyoff + izoff;
              CONTRACT_VAR(ncells);
              CHECK(icell < ncells);
              values[icell] += value*W(etaMag, Hdeti)*Hdeti;
            }
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// The straightforward binning method.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
vector<Value>
binFieldList2Lattice(const FieldList<Dimension, Value>& fieldList,
                     const typename Dimension::Vector& xmin,
                     const typename Dimension::Vector& xmax,
                     const vector<unsigned>& nsample) {

  // Preconditions.
  REQUIRE(nsample.size() == Dimension::nDim);

  // Some convenient typedefs.
  typedef typename Dimension::Vector Vector;

  // Get the parallel geometry.
  //const unsigned procID = Process::getRank();
  const unsigned numProcs = Process::getTotalNumberOfProcesses();
  CONTRACT_VAR(numProcs);
  CHECK(numProcs >= 1);

  // We need to exclude any nodes that come from the Distributed boundary condition.
#ifdef SPHERAL_ENABLE_MPI
  BoundingVolumeDistributedBoundary<Dimension>& distributedBoundary = BoundingVolumeDistributedBoundary<Dimension>::instance();
#endif

  // Compute the total number of sample points.
  unsigned ntotal = 1;
  for (size_t i = 0; i != Dimension::nDim; ++i) ntotal *= nsample[i];
  CHECK(ntotal > 0);

  // Prepare the result.
  vector<Value> result(ntotal, DataTypeTraits<Value>::zero());

  // Loop over the nodes.
  for (AllNodeIterator<Dimension> nodeItr = fieldList.nodeBegin();
       nodeItr != fieldList.nodeEnd();
       ++nodeItr) {

    // Mask out any parallel boundary nodes.
#ifdef SPHERAL_ENABLE_MPI
    const bool useNode = ((not distributedBoundary.haveNodeList(*nodeItr.nodeListPtr())) or
                          count(distributedBoundary.ghostNodes(*(nodeItr.nodeListPtr())).begin(),
                                distributedBoundary.ghostNodes(*(nodeItr.nodeListPtr())).end(),
                                nodeItr.nodeID()) == 0);
#else
    const bool useNode = true;
#endif
    const Vector& xi = nodeItr.position();
    if (useNode and testPointInBox(xi, xmin, xmax)) {

      // Add this nodes contribution to it's cell.
      const size_t i = latticeIndex(xi, xmin, xmax, nsample);
      result[i] += fieldList(nodeItr);

    }
  }

  // In parallel we have to reduce the elements across processors.
#ifdef SPHERAL_ENABLE_MPI
  vector<char> localBuffer;
  packElement(result, localBuffer);
  int bufSize = localBuffer.size();
  
  // This is done so that everyone sums stuff up in the same order and we should
  // get the same answer on each processor to bit precision.
  result = vector<Value>(ntotal, DataTypeTraits<Value>::zero());

  // Check that everyone agrees about the size.
  CHECK(bufSize == allReduce(bufSize, SPHERAL_OP_MIN));
  CHECK(bufSize == allReduce(bufSize, SPHERAL_OP_MAX));

  // Sum up everyone's contribution.
  for (auto sendProc = 0u; sendProc != numProcs; ++sendProc) {
    vector<char> buffer(localBuffer);
    MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
    vector<Value> localResult;
    vector<char>::const_iterator itr = buffer.begin();
    unpackElement(localResult, itr, buffer.end());
    CHECK(itr == buffer.end());
    CHECK(localResult.size() == result.size());
    for (size_t i = 0; i != ntotal; ++i) result[i] += localResult[i];
  }
#endif

  // We're done.
  return result;
}

//------------------------------------------------------------------------------
// Bin with a smoothing applied.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
vector<Value>
binFieldList2Lattice(const FieldList<Dimension, Value>& fieldList,
                     const TableKernel<Dimension>& W,
                     const typename Dimension::Vector& xmin,
                     const typename Dimension::Vector& xmax,
                     const vector<unsigned>& nsample) {

  // Preconditions.
  REQUIRE(nsample.size() == Dimension::nDim);

  // Get the parallel geometry.
  //const unsigned procID = Process::getRank();
  const unsigned numProcs = Process::getTotalNumberOfProcesses();
  CONTRACT_VAR(numProcs);
  CHECK(numProcs >= 1);

  // We need to exclude any nodes that come from the Distributed boundary condition.
#ifdef SPHERAL_ENABLE_MPI
  BoundingVolumeDistributedBoundary<Dimension>& distributedBoundary = BoundingVolumeDistributedBoundary<Dimension>::instance();
#endif

  // Compute the total number of sample points.
  unsigned ntotal = 1;
  for (size_t i = 0; i != Dimension::nDim; ++i) ntotal *= nsample[i];
  CHECK(ntotal > 0);

  // Prepare the result.
  vector<Value> result(ntotal, DataTypeTraits<Value>::zero());

  // Loop over the nodes.
  for (AllNodeIterator<Dimension> nodeItr = fieldList.nodeBegin();
       nodeItr != fieldList.nodeEnd();
       ++nodeItr) {

    // Mask out any parallel boundary nodes.
#ifdef SPHERAL_ENABLE_MPI
    const bool useNode = ((not distributedBoundary.haveNodeList(*nodeItr.nodeListPtr())) or
                          count(distributedBoundary.ghostNodes(*(nodeItr.nodeListPtr())).begin(),
                                distributedBoundary.ghostNodes(*(nodeItr.nodeListPtr())).end(),
                                nodeItr.nodeID()) == 0);
#else
    const bool useNode = true;
#endif
    if (useNode and testPointInBox(nodeItr.position(), xmin, xmax)) {

      // Add this nodes contribution to it's cells.
      incrementCellValues(result, fieldList(nodeItr), nodeItr.position(), nodeItr.H(), xmin, xmax, nsample, W);

    }
  }

  // In parallel we have to reduce the elements across processors.
#ifdef SPHERAL_ENABLE_MPI
  vector<char> localBuffer;
  packElement(result, localBuffer);
  int bufSize = localBuffer.size();
  
  // This is done so that everyone sums stuff up in the same order and we should
  // get the same answer on each processor to bit precision.
  result = vector<Value>(ntotal, DataTypeTraits<Value>::zero());

  BEGIN_CONTRACT_SCOPE
  // Check that everyone agrees about the size.
  {
    CHECK(bufSize == allReduce(bufSize, SPHERAL_OP_MIN));
  }
  END_CONTRACT_SCOPE

  // Sum up everyone's contribution.
  for (auto sendProc = 0u; sendProc != numProcs; ++sendProc) {
    vector<char> buffer(localBuffer);
    MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
    vector<Value> localResult;
    vector<char>::const_iterator itr = buffer.begin();
    unpackElement(localResult, itr, buffer.end());
    CHECK(itr == buffer.end());
    CHECK(localResult.size() == result.size());
    for (size_t i = 0; i != ntotal; ++i) result[i] += localResult[i];
  }
#endif

  // We're done.
  return result;
}

}
