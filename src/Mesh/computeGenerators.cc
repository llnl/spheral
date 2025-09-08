//---------------------------------Spheral++----------------------------------//
// computeGenerators.
//
// Helper method used when we're generating meshes.
// This method takes a set of NodeLists, figures out what domains have to talk
// to one another, and returns the flattened set of positions and Hs for the 
// the generators this domain needs (including those from neighbor domains).
//
// Created by JMO, Mon Dec  6 10:34:43 PST 2010
//----------------------------------------------------------------------------//
#include "computeGenerators.hh"
#include "Mesh.hh"
#include "MeshConstructionUtilities.hh"
#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "Neighbor/Neighbor.hh"
#include "Boundary/Boundary.hh"
#include "Utilities/bisectSearch.hh"
#include "Utilities/packElement.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/removeElements.hh"
#include "Distributed/Communicator.hh"

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <algorithm>
#include <set>
using std::vector;
using std::set;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// The method itself.
//------------------------------------------------------------------------------
template<typename Dimension, typename NodeListIterator, typename BoundaryIterator>
void
computeGenerators(NodeListIterator nodeListBegin,
                  NodeListIterator nodeListEnd,
                  BoundaryIterator boundaryBegin,
                  BoundaryIterator boundaryEnd,
                  const bool meshGhostNodes,
                  const typename Dimension::Vector& xmin,
                  const typename Dimension::Vector& xmax,
                  vector<typename Dimension::Vector>& positions,
                  vector<typename Dimension::SymTensor>& Hs,
                  vector<unsigned>& offsets) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;

  // Parallel geometry.
  const unsigned numNodeLists = distance(nodeListBegin, nodeListEnd);

  // Flatten the local positions and Hs.
  vector<Vector> localPositions;
  vector<SymTensor> localHs;
  offsets = vector<unsigned>(1, 0);
  unsigned i, k, nlocal = 0;
  for (NodeListIterator nodeListItr = nodeListBegin; nodeListItr != nodeListEnd; ++nodeListItr) {

    // We always want all the internal nodes.
    const Field<Dimension, Vector>& pos = (**nodeListItr).positions();
    const Field<Dimension, SymTensor>& H = (**nodeListItr).Hfield();
    copy(pos.internalBegin(), pos.internalEnd(), back_inserter(localPositions));
    copy(H.internalBegin(), H.internalEnd(), back_inserter(localHs));
    offsets.push_back(offsets.back() + (**nodeListItr).numInternalNodes());
    nlocal += (**nodeListItr).numInternalNodes();

    // Look for any boundaries whose ghost nodes we should include.
    if (meshGhostNodes) {
      for (BoundaryIterator bcItr = boundaryBegin;
           bcItr != boundaryEnd;
           ++bcItr) {
        if ((*bcItr)->meshGhostNodes() and (*bcItr)->haveNodeList(**nodeListItr)) {
          const auto& ghostNodes = (*bcItr)->ghostNodes(**nodeListItr);
          for (k = 0; k < ghostNodes.size(); ++k) {
            i = ghostNodes[k];
            //           if (!(i < pos.numElements())) {
            //             cerr << "Blago!  " << ghostNodes.size() << endl << "  ---> ";
            //             for (unsigned kk = 0; kk != ghostNodes.size(); ++kk) cerr << ghostNodes[kk] << " ";
            //             cerr << endl;
            //           }
            CHECK2(i < pos.numElements(), "Out of bounds:  " << i << " " << pos.numElements());
            if (testPointInBox(pos(i), xmin, xmax)) {
              ++offsets.back();
              ++nlocal;
              localPositions.push_back(pos(i));
              localHs.push_back(H(i));
            }
          }
        }
      }
    }
  }

  // // Look for any boundaries whose ghost nodes we should include.
  // if (meshGhostNodes) {
  //   for (BoundaryIterator bcItr = boundaryBegin;
  //        bcItr != boundaryEnd;
  //        ++bcItr) {
  //     for (NodeListIterator nodeListItr = nodeListBegin; nodeListItr != nodeListEnd; ++nodeListItr) {
  //       if ((*bcItr)->meshGhostNodes() and (*bcItr)->haveNodeList(**nodeListItr)) {
  //         const Field<Dimension, Vector>& pos = (**nodeListItr).positions();
  //         const Field<Dimension, SymTensor>& H = (**nodeListItr).Hfield();
  //         const vector<int>& ghostNodes = (*bcItr)->ghostNodes(**nodeListItr);
  //         for (k = 0; k != ghostNodes.size(); ++k) {
  //           i = ghostNodes[k];
  //           //           if (!(i < pos.numElements())) {
  //           //             cerr << "Blago!  " << ghostNodes.size() << endl << "  ---> ";
  //           //             for (unsigned kk = 0; kk != ghostNodes.size(); ++kk) cerr << ghostNodes[kk] << " ";
  //           //             cerr << endl;
  //           //           }
  //           CHECK2(i < pos.numElements(), "Out of bounds:  " << i << " " << pos.numElements());
  //           if (testPointInBox(pos(i), xmin, xmax)) {
  //             ++offsets.back();
  //             ++nlocal;
  //             localPositions.push_back(pos(i));
  //             localHs.push_back(H(i));
  //           }
  //         }
  //       }
  //     }
  //   }
  // }

  CONTRACT_VAR(numNodeLists);
  CHECK(localPositions.size() == nlocal);
  CHECK(localHs.size() == nlocal);
  CHECK(offsets.size() == numNodeLists + 1);

  // Copy our local info to the result arrays.
  positions = localPositions;
  Hs = localHs;

  // That's it.
}

}

