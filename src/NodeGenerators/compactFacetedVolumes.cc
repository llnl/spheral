#include <algorithm>
#include <ctime>

#include "compactFacetedVolumes.hh"

#include "NodeList/FluidNodeList.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/GammaLawGas.hh"
#include "DataBase/DataBase.hh"
#include "Utilities/packElement.hh"
#include "Distributed/allReduce.hh"

#include "config.hh"

#ifdef ENABLE_MPI
#include <mpi.h>
#endif

namespace Spheral {


//------------------------------------------------------------------------------
// Push FacetedVolume shapes together inside a surface, but excluding mutual
// overlap.
//------------------------------------------------------------------------------
template<typename Dimension>
unsigned compactFacetedVolumes(std::vector<typename Dimension::FacetedVolume>& shapes,
                               std::vector<typename Dimension::Vector>& centers,
                               std::vector<int>& flags,
                               const typename Dimension::FacetedVolume& surface,
                               const double depthmax,
                               const unsigned surfaceIterations,
                               const unsigned maxIterations,
                               const double dispfrac,
                               const double maxoverlapfrac) {

  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  // Preconditions.
  const unsigned nshapes = shapes.size();
  VERIFY(centers.size() == nshapes);
  VERIFY(flags.size() == nshapes);

  // Only proceed if there's work to do!
  int flagmax = *max_element(flags.begin(), flags.end());
  if (allReduce(flagmax, SPHERAL_OP_MAX) != 2) return 0;

  // Carve up the shapes range in parallel.
  // const size_t ndomain0 = nshapes/nprocs;
  // const size_t remainder = nshapes % nprocs;
  // CHECK(remainder < nprocs);
  // size_t ndomain = ndomain0;
  // if (rank < remainder) ndomain += 1;
  // const size_t imin = rank*ndomain0 + min(rank, remainder);
  // const size_t imax = imin + ndomain;
  size_t imin, imax;
  if (Process::getRank() == 0) {
    imin = 0;
    imax = nshapes;
  } else {
    imin = 0;
    imax = 0;
  }

  // Build a fake NodeList with Neighbor to help selecting shapes that might be neighbors.
  const double length = (surface.xmax() - surface.xmin()).maxElement();
  const Vector xghost = surface.xmax() + length*Vector::one;

  // Make a temporary NodeList so we can use it's neighbor logic.
  PhysicalConstants constants(1.0, 1.0, 1.0);
  GammaLawGas<Dimension> eos(2.0, 2.0, constants, 0.0, 1e10, MaterialPressureMinType::PressureFloor, 0.0);
  FluidNodeList<Dimension> nodes("shapes", eos, nshapes, 0,
                                            1e-10, 1e10, 1.0, 1.01, 500, 0.0, 1e10);
  NestedGridNeighbor<Dimension> neighbor(nodes, NeighborSearchType::GatherScatter, 31, 3.0*length, surface.xmin(), 1.0, 1);
  nodes.registerNeighbor(neighbor);
  DataBase<Dimension> db;
  db.appendNodeList(nodes);
  Field<Dimension, Vector>& pos = nodes.positions();
  Field<Dimension, SymTensor>& H = nodes.Hfield();
  Field<Dimension, double> radius("radius", nodes);

  // Scratch fields for Neighbor operations.
  vector<int> masterList, coarseNeighbors, refineNeighbors;

  // Figure out the effective size per shape.
  for (size_t i = 0; i < nshapes; ++i) {
    if (flags[i] != 0) {
      pos[i] = centers[i];
    } else {
      pos[i] = xghost;
    }
    double Ri = Dimension::rootnu(shapes[i].volume()/M_PI);
    radius[i] = Ri;
    H[i] = SymTensor::one / (2.0*Ri);
  }

  // First iterate the requested number of iterations with surface repulsion.
  unsigned iter = 0;
  double maxoverlap = 10.0*maxoverlapfrac;
  std::clock_t
    t0,
    tsurface = std::clock_t(0),
    tothers = std::clock_t(0),
    toverlap = std::clock_t(0);


  if (Process::getRank() == 0) {  // Fastest in serial
    while (iter < surfaceIterations or (iter < maxIterations and maxoverlap > maxoverlapfrac)) {
      iter += 1;
      vector<Vector> displacements(nshapes);

      // Update the neighbor info.
      neighbor.updateNodes();

      // Add the repulsive component from the asteroid surface.
      t0 = std::clock();
      for (auto i = imin; i < imax; ++i) {
        if (flags[i] == 2) {
          const Vector p = surface.closestPoint(centers[i]);
          const double d = (p - centers[i]).magnitude(); //  surface.distance(centers[i]);
          if (surface.contains(centers[i])) {
            displacements[i] += dispfrac*(1.0 - double(iter)/surfaceIterations)*max(0.01, 1.0 - d/depthmax)*Dimension::rootnu(shapes[i].volume()/M_PI)*(centers[i] - p).unitVector();
          } else {
            displacements[i] -= dispfrac*(1.0 - double(iter)/surfaceIterations)*Dimension::rootnu(shapes[i].volume()/M_PI)*(centers[i] - p).unitVector();
          }
        }
      }
      tsurface += std::clock() - t0;

      // Look for any overlapping shapes and drive them apart.
      t0 = std::clock();
      for (auto i = imin; i < imax; ++i) {
        if (flags[i] == 3) {
          const FacetedVolume bi = shapes[i] + centers[i];
          neighbor.setMasterList(i, masterList, coarseNeighbors);
          neighbor.setRefineNeighborList(i, coarseNeighbors, refineNeighbors);
          for (auto jitr = refineNeighbors.begin(); jitr != refineNeighbors.end(); ++jitr) {
            const auto j = *jitr;
            if (j != (int)i) {
              if ((flags[j] == 1 and bi.intersect(shapes[j])) or
                  (j > (int)i and flags[j] >= 2 and bi.intersect(shapes[j] + centers[j]))) {
                const double Vi = shapes[i].volume();
                const double Vj = shapes[j].volume();
                const double ri = Dimension::rootnu(Vi/M_PI);
                const double rj = Dimension::rootnu(Vj/M_PI);
                const Vector dji = centers[j] - centers[i];
                const double djimag = dji.magnitude();
                const Vector delta = dji.unitVector() * dispfrac*max(0.1*min(ri, rj), max(0.0, ri + rj - djimag));
                if (flags[j] == 3) {
                  displacements[i] -= 0.25*Vi/(Vi + Vj)*delta;
                  displacements[j] += 0.25*Vj/(Vi + Vj)*delta;
                } else {
                  displacements[i] -= 0.25*delta;
                }
              }
            }
          }
        }
      }

      // Apply the displacements of this iteration.
      for (auto i = imin; i < imax; ++i) {
        if (flags[i] >= 2) {
          centers[i] += displacements[i];
          pos[i] = centers[i];
        }
      }
      tothers += std::clock() - t0;

      // #ifdef ENABLE_MPI
      //     // Global broadcast of the new centers.
      //     for (int iproc = 0; iproc != nprocs; ++iproc) {
      //       int jmin = imin, jmax = imax;
      //       MPI_Bcast(&jmin, 1, MPI_INT, iproc, Communicator::communicator());
      //       MPI_Bcast(&jmax, 1, MPI_INT, iproc, Communicator::communicator());
      //       vector<char> buffer;
      //       int bufsize = 0;
      //       if (rank == iproc) {
      //         for (auto i = imin; i != imax; ++i) packElement(centers[i], buffer);
      //         bufsize = buffer.size();
      //       }
      //       MPI_Bcast(&bufsize, 1, MPI_INT, iproc, Communicator::communicator());
      //       buffer.resize(bufsize);
      //       MPI_Bcast(&buffer[0], bufsize, MPI_CHAR, iproc, Communicator::communicator());
      //       if (rank != iproc) {
      //         vector<char>::const_iterator bufitr = buffer.begin();
      //         for (auto j = jmin; j < jmax; ++j) unpackElement(centers[j], bufitr, buffer.end());
      //         CHECK(bufitr == buffer.end());
      //       }
      //     }
      // #endif

      //  Check the current level of overlap.
      t0 = std::clock();
      neighbor.updateNodes();
      maxoverlap = 0.0;
      for (auto i = imin; i < imax; ++i) {
        if (flags[i] >= 2) {
          flags[i] = 2;
          const FacetedVolume shapei = shapes[i] + centers[i];
          const double Ri = radius[i];
          // const vector<vector<int>>& neighbors = cm.connectivityForNode(0, i);
          // CHECK(neighbors.size() == 1);
          // for (auto j: neighbors[0]) {
          // for (auto j = 0; j != nshapes; ++j) {
          neighbor.setMasterList(i, masterList, coarseNeighbors);
          neighbor.setRefineNeighborList(i, coarseNeighbors, refineNeighbors);
          for (auto jitr = refineNeighbors.begin(); jitr != refineNeighbors.end(); ++jitr) {
            const auto j = *jitr;
            if (j != (int)i) {
              if (flags[j] >= 1) {
                FacetedVolume shapej;
                if (flags[j] == 1) {
                  shapej = shapes[j];
                } else {
                  shapej = shapes[j] + centers[j];
                }
                if (shapei.intersect(shapej)) {
                  const double Rj = radius[j];
                  auto overlap = max(0.0, (Ri + Rj - (shapei.centroid() - shapej.centroid()).magnitude())/(Ri + Rj));

                  // # Ri = (shapei.volume/pi)**(1.0/self.ndim)
                  // # Rj = (shapej.volume/pi)**(1.0/self.ndim)
                  // # vol = min(shapei.volume, shapej.volume)
                  // # verticesi = shapei.vertices()
                  // # verticesj = shapej.vertices()
                  // # dmax = 0.0
                  // # for vi in verticesi:
                  // #     if shapej.contains(vi):
                  // #         dmax = max(dmax, shapej.distance(vi))
                  // # for vj in verticesj:
                  // #     if shapei.contains(vj):
                  // #         dmax = max(dmax, shapei.distance(vj))
                  // # overlap = dmax/(vol/pi)**(1.0/self.ndim)

                  maxoverlap = max(maxoverlap, overlap);
                  if (overlap > maxoverlapfrac) flags[i] = 3;
                }
              }
            }
          }
          if (flags[i] != 3 and iter >= surfaceIterations) {
            // We can freeze this shape now.
            flags[i] = 1;
            shapes[i] += centers[i];
            pos[i] = centers[i];
          }
        }
      }
      toverlap += std::clock() - t0;

      // #ifdef ENABLE_MPI
      //     // Global broadcast of the new flags.
      //     for (int iproc = 0; iproc != nprocs; ++iproc) {
      //       int jmin = imin, jmax = imax;
      //       MPI_Bcast(&jmin, 1, MPI_INT, iproc, Communicator::communicator());
      //       MPI_Bcast(&jmax, 1, MPI_INT, iproc, Communicator::communicator());
      //       vector<char> buffer;
      //       int bufsize = 0;
      //       if (rank == iproc) {
      //         for (auto i = imin; i != imax; ++i) {
      //           packElement(flags[i], buffer);
      //           packElement(centers[i], buffer);
      //           packElement(shapes[i], buffer);
      //         }
      //         bufsize = buffer.size();
      //       }
      //       MPI_Bcast(&bufsize, 1, MPI_INT, iproc, Communicator::communicator());
      //       buffer.resize(bufsize);
      //       MPI_Bcast(&buffer[0], bufsize, MPI_CHAR, iproc, Communicator::communicator());
      //       if (rank != iproc) {
      //         vector<char>::const_iterator bufitr = buffer.begin();
      //         for (auto j = jmin; j < jmax; ++j) {
      //           unpackElement(flags[j], bufitr, buffer.end());
      //           unpackElement(centers[j], bufitr, buffer.end());
      //           unpackElement(shapes[j], bufitr, buffer.end());
      //         }
      //         CHECK(bufitr == buffer.end());
      //       }
      //     }
      //     maxoverlap = allReduce(maxoverlap, SPHERAL_OP_MAX);
      // #endif
      // double sumdisp = allReduce(std::accumulate(displacements.begin(), displacements.end(), 0.0, [](const double prior, const Vector& elemval) { return prior + elemval.magnitude(); }), SPHERAL_OP_SUM);
      // if (rank == 0) {
      //   cout << "   Iteration " << iter 
      //        << ", maxoverlap " << maxoverlap 
      //        << ", total 1's " << std::accumulate(flags.begin(), flags.end(), 0, [](const int prior, const int elemval) { return elemval == 1 ? prior + 1 : prior; })
      //        << ", total 2's " << std::accumulate(flags.begin(), flags.end(), 0, [](const int prior, const int elemval) { return elemval == 2 ? prior + 1 : prior; })
      //        << ", total 3's " << std::accumulate(flags.begin(), flags.end(), 0, [](const int prior, const int elemval) { return elemval == 3 ? prior + 1 : prior; })
      //        << ", sum displacements " << sumdisp
      //        << ", fraction of stopping criteria " << maxoverlapfrac/maxoverlap << endl;
      // }

    } // end of iteration
  }
  iter = allReduce(iter, SPHERAL_OP_MAX);

  // Any shapes we were unable to disentangle turn back to inactive, otherwise set the successful
  // survivors to flag=1.
  for (auto i = imin; i < imax; ++i) {
    if (flags[i] == 3) {
      CHECK(maxoverlap > maxoverlapfrac);
      flags[i] = 0;
      centers[i] = Vector::zero;
      // # # Make one last ditch attempt to randomly fit this shape in.
      // # centers[i] = self.randomCenter(i, centers)
      // # if centers[i]:
      // #     flags[i] = 1
      // # else:
      // #     flags[i] = 0
    }
  }

#ifdef ENABLE_MPI
  const size_t rank = Process::getRank();
  const size_t nprocs = Process::getTotalNumberOfProcesses();
  // Global broadcast of the new geometry.
  for (auto iproc = 0u; iproc != nprocs; ++iproc) {
    int jmin = imin, jmax = imax;
    MPI_Bcast(&jmin, 1, MPI_INT, iproc, Communicator::communicator());
    MPI_Bcast(&jmax, 1, MPI_INT, iproc, Communicator::communicator());
    vector<char> buffer;
    int bufsize = 0;
    if (rank == iproc) {
      for (auto i = imin; i != imax; ++i) {
        packElement(flags[i], buffer);
        packElement(centers[i], buffer);
        packElement(shapes[i], buffer);
      }
      bufsize = buffer.size();
    }
    MPI_Bcast(&bufsize, 1, MPI_INT, iproc, Communicator::communicator());
    buffer.resize(bufsize);
    MPI_Bcast(&buffer[0], bufsize, MPI_CHAR, iproc, Communicator::communicator());
    if (rank != iproc) {
      vector<char>::const_iterator bufitr = buffer.begin();
      for (auto j = jmin; j < jmax; ++j) {
        unpackElement(flags[j], bufitr, buffer.end());
        unpackElement(centers[j], bufitr, buffer.end());
        unpackElement(shapes[j], bufitr, buffer.end());
      }
      CHECK(bufitr == buffer.end());
    }
  }
#endif

  // That's it.
  if (Process::getRank() == 0) cout << "compactFacetedVolume timing:" 
                                    << " tsurface=" << (tsurface / (double) CLOCKS_PER_SEC) 
                                    << " tothers=" << (tothers / (double) CLOCKS_PER_SEC) 
                                    << " toverlap=" << (toverlap / (double) CLOCKS_PER_SEC)  << endl;
  return iter;
}

}
