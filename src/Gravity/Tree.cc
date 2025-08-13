//---------------------------------Spheral++----------------------------------//
// Tree -- a implementation of a quad/oct (2D/3D) tree structure.
//
// Extracted from our original implementation in TreeGravity.
// This class carries along a few data members from that heritage (like the mass
// and velocity) which we retain for convenience. Those attributes can be ignored
// for purely geometrical applications.
//
// Created by JMO, Tue Oct  4 10:17:41 PDT 2022
//----------------------------------------------------------------------------//
#include "Tree.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/packElement.hh"
#include "Distributed/allReduce.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/PairComparisons.hh"
#include "Distributed/Communicator.hh"
#include "Utilities/DBC.hh"
#include "Utilities/Hashes.hh"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>
using std::vector;
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
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
Tree<Dimension>::
Tree(const Vector& xmin, const Vector& xmax):
  mBoxLength((xmax - xmin).maxAbsElement()),
  mXmin(xmin),
  mXmax(xmax),
  mLevels() {
  REQUIRE(mBoxLength > 0.0);
  REQUIRE(mXmin.x() <= mXmax.x());
  REQUIRE(mXmin.y() <= mXmax.y());
  REQUIRE(mXmin.z() <= mXmax.z());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Tree<Dimension>::
~Tree() {
}

//------------------------------------------------------------------------------
// dumpTree
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
Tree<Dimension>::
dumpTree(const bool globalTree) const {
  std::stringstream ss;
  CellKey ix, iy, iz;
#ifdef USE_MPI
  const unsigned numProcs = Process::getTotalNumberOfProcesses();
  const unsigned rank = Process::getRank();
#endif
  unsigned nlevels = mLevels.size();
  if (globalTree) nlevels = allReduce(nlevels, SPHERAL_OP_MAX);

  ss << "Tree : nlevels = " << nlevels << "\n";
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {

    // Gather up the level cells and sort them.
    vector<Cell> cells;
    vector<char> localBuffer;
    cells.reserve(mLevels[ilevel].size());
    if (ilevel < mLevels.size()) {
      for (typename TreeLevel::const_iterator itr = mLevels[ilevel].begin();
           itr != mLevels[ilevel].end();
           ++itr) {
        cells.push_back(itr->second);
        this->serialize(itr->second, localBuffer);
      }
    }
#ifdef USE_MPI
    if (globalTree) {
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
        if (bufSize > 0) {
          vector<char> buffer = localBuffer;
          MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
          if (rank != sendProc) {
            vector<char>::const_iterator itr = buffer.begin();
            while (itr < buffer.end()) {
              cells.push_back(Cell());
              this->deserialize(cells.back(), itr, buffer.end());
            }
          }
        }
      }
    }
#endif
    sort(cells.begin(), cells.end());
    cells.erase(unique(cells.begin(), cells.end()), cells.end());

    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : numCells = " << cells.size() << "\n";
    for (typename vector<Cell>::const_iterator itr = cells.begin();
         itr != cells.end();
         ++itr) {
      const Cell& cell = *itr;
      extractCellIndices(cell.key, ix, iy, iz);
      ss << "    Cell key=" << cell.key << " : (ix,iy,iz)=(" << ix << " " << iy << " " << iz << ")\n"
         << "         xcm=" << cell.xcm << " vcm=" << cell.vcm << " rcm2cc=" << sqrt(cell.rcm2cc2) << " M=" << cell.M  << " Mglobal=" << cell.Mglobal << "\n"
         << "         daughters = ( ";
      for (vector<CellKey>::const_iterator ditr = cell.daughters.begin();
           ditr != cell.daughters.end();
           ++ditr) ss << *ditr << " ";
      ss << ")\n"
         << "         nodes = [";
      for (unsigned k = 0; k != cell.masses.size(); ++k) ss << " ("
                                                            << cell.masses[k] << " "
                                                            << cell.positions[k] << " "
                                                            << cell.velocities[k] << ")";
      ss <<" ]\n";
    }
  }
  return ss.str();
}

//------------------------------------------------------------------------------
// dumpTreeStatistics
//------------------------------------------------------------------------------
template<typename Dimension>
std::string
Tree<Dimension>::
dumpTreeStatistics(const bool globalTree) const {
  std::stringstream ss;
#ifdef USE_MPI
  const unsigned numProcs = Process::getTotalNumberOfProcesses();
  const unsigned rank = Process::getRank();
#endif
  unsigned nlevels = mLevels.size();
  if (globalTree) nlevels = allReduce(nlevels, SPHERAL_OP_MAX);

  ss << "Tree : nlevels = " << nlevels << "\n";
  for (unsigned ilevel = 0; ilevel != nlevels; ++ilevel) {

    // Gather up the level cells and sort them.
    vector<Cell> cells;
    vector<char> localBuffer;
    cells.reserve(mLevels[ilevel].size());
    if (ilevel < mLevels.size()) {
      for (typename TreeLevel::const_iterator itr = mLevels[ilevel].begin();
           itr != mLevels[ilevel].end();
           ++itr) {
        cells.push_back(itr->second);
        this->serialize(itr->second, localBuffer);
      }
    }
#ifdef USE_MPI
    if (globalTree) {
      for (unsigned sendProc = 0; sendProc != numProcs; ++sendProc) {
        unsigned bufSize = localBuffer.size();
        MPI_Bcast(&bufSize, 1, MPI_UNSIGNED, sendProc, Communicator::communicator());
        if (bufSize > 0) {
          vector<char> buffer = localBuffer;
          MPI_Bcast(&buffer.front(), bufSize, MPI_CHAR, sendProc, Communicator::communicator());
          if (rank != sendProc) {
            vector<char>::const_iterator itr = buffer.begin();
            while (itr < buffer.end()) {
              cells.push_back(Cell());
              this->deserialize(cells.back(), itr, buffer.end());
            }
          }
        }
      }
    }
#endif
    sort(cells.begin(), cells.end());
    cells.erase(unique(cells.begin(), cells.end()), cells.end());

    // Count up how much of everything we have.
    ss << "--------------------------------------------------------------------------------\n" 
       << " Level " << ilevel << " : numCells = " << cells.size() << "\n";
    unsigned nparticles = 0;
    for (typename vector<Cell>::const_iterator itr = cells.begin();
         itr != cells.end();
         ++itr) {
      const Cell& cell = *itr;
      nparticles += cell.masses.size();
    }
    ss << "         : nparts = " << nparticles << "\n";
  }
  return ss.str();
}

//------------------------------------------------------------------------------
// Serialize a tree to a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Tree<Dimension>::
serialize(std::vector<char>& buffer) const {
  const auto nlevels = this->numLevels();
  packElement(nlevels, buffer);
  for (auto ilevel = 0u; ilevel < nlevels; ++ilevel) {
    const auto& cells = mLevels[ilevel];
    const auto ncells = cells.size();
    packElement(ncells, buffer);
    for (typename TreeLevel::const_iterator itr = cells.begin();
         itr != cells.end();
         ++itr) {
      packElement(itr->first, buffer);
      serialize(itr->second, buffer);
    }
  }
}

//------------------------------------------------------------------------------
// Serialize a cell to a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Tree<Dimension>::
serialize(const Tree<Dimension>::Cell& cell,
          std::vector<char>& buffer) const {
  packElement(cell.M, buffer);
  packElement(cell.Mglobal, buffer);
  packElement(cell.xcm, buffer);
  packElement(cell.vcm, buffer);
  packElement(cell.rcm2cc2, buffer);
  packElement(cell.key, buffer);
  packElement(cell.daughters, buffer);
  packElement(cell.masses, buffer);
  packElement(cell.positions, buffer);
  packElement(cell.velocities, buffer);
}

//------------------------------------------------------------------------------
// Deserialize a tree from a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Tree<Dimension>::
deserialize(vector<char>::const_iterator& bufItr,
            const vector<char>::const_iterator& endItr) {
  size_t nlevels;
  CellKey key;
  Cell cell;
  unpackElement(nlevels, bufItr, endItr);
  this->numLevels(nlevels);
  for (auto ilevel = 0u; ilevel < nlevels; ++ilevel) {
    size_t ncells;
    unpackElement(ncells, bufItr, endItr);
    for (auto i = 0u; i < ncells; ++i) {
      cell.daughters = vector<CellKey>();
      cell.masses = vector<double>();
      cell.positions = vector<Vector>();
      cell.velocities = vector<Vector>();
      unpackElement(key, bufItr, endItr);
      deserialize(cell, bufItr, endItr);
      mLevels[ilevel][key] = cell;
    }
  }
  constructDaughterPtrs();
}

//------------------------------------------------------------------------------
// Deserialize a cell from a buffer of char.
//------------------------------------------------------------------------------
template<typename Dimension>
void
Tree<Dimension>::
deserialize(Tree<Dimension>::Cell& cell,
            vector<char>::const_iterator& bufItr,
            const vector<char>::const_iterator& endItr) const {
  unpackElement(cell.M, bufItr, endItr);
  unpackElement(cell.Mglobal, bufItr, endItr);
  unpackElement(cell.xcm, bufItr, endItr);
  unpackElement(cell.vcm, bufItr, endItr);
  unpackElement(cell.rcm2cc2, bufItr, endItr);
  unpackElement(cell.key, bufItr, endItr);
  unpackElement(cell.daughters, bufItr, endItr);
  unpackElement(cell.masses, bufItr, endItr);
  unpackElement(cell.positions, bufItr, endItr);
  unpackElement(cell.velocities, bufItr, endItr);
}

}
