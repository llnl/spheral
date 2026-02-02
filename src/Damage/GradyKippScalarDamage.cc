//---------------------------------Spheral++----------------------------------//
// GradyKippScalarDamage -- Implements a scalar damage model on 
// SolidNodeLists based on the Grady & Kipp crack model.
//
// References:
//   Benz, W. & Asphaug, E., 1995 "Computer Physics Comm.", 87, 253-265.
//   Benz, W. & Asphaug, E., 1994 "Icarus", 107, 98-116.
//
// Created by JMO, Mon Sep 20 21:48:53 PDT 2004
//----------------------------------------------------------------------------//
#include <string>

#include "GradyKippScalarDamage.hh"
#include "SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "ScalarDamagePolicy.hh"
#include "weibullFlawDistribution.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Kernel/TableKernel.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FileIO/FileIO.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GradyKippScalarDamage<Dimension>::
GradyKippScalarDamage(SolidNodeList<Dimension>& nodeList,
                      FluidNodeList<Dimension>& damagedNodeList,
                      const double kWeibull,
                      const double mWeibull,
                      const double volume,
                      const TableKernel<Dimension>& kernel,
                      const unsigned seed,
                      const double crackGrowthMultiplier,
                      const int minFlawsPerNode,
                      const int minTotalFlaws,
                      const double averageFailure):
  ScalarDamageModel<Dimension>(nodeList,
                               damagedNodeList,
                               kernel.kernelExtent(),
                               crackGrowthMultiplier,
                               weibullFlawDistribution(volume,
                                                       seed,
                                                       kWeibull,
                                                       mWeibull,
                                                       nodeList,
                                                       minFlawsPerNode,
                                                       minTotalFlaws,
                                                       averageFailure)),
  mkWeibull(kWeibull),
  mmWeibull(mWeibull),
  mSeed(seed) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GradyKippScalarDamage<Dimension>::
~GradyKippScalarDamage() {
}

//------------------------------------------------------------------------------
// Access the internal parameters of the model.
//------------------------------------------------------------------------------
template<typename Dimension>
double
GradyKippScalarDamage<Dimension>::
kWeibull() const {
  return mkWeibull;
}

template<typename Dimension>
double
GradyKippScalarDamage<Dimension>::
mWeibull() const {
  return mmWeibull;
}

template<typename Dimension>
unsigned
GradyKippScalarDamage<Dimension>::
seed() const {
  return mSeed;
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GradyKippScalarDamage<Dimension>::
dumpState(FileIO& file, const string& pathName) const {

  ScalarDamageModel<Dimension>::dumpState(file, pathName);

  file.write(mkWeibull, pathName + "/kWeibull");
  file.write(mmWeibull, pathName + "/mWeibull");
  file.write((int) mSeed, pathName + "/seed");

}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GradyKippScalarDamage<Dimension>::
restoreState(const FileIO& file, const string& pathName) {

  ScalarDamageModel<Dimension>::restoreState(file, pathName);

  file.read(mkWeibull, pathName + "/kWeibull");
  file.read(mmWeibull, pathName + "/mWeibull");
  int seed;
  file.read(seed, pathName + "/seed");
  mSeed = seed;
}

}

//------------------------------------------------------------------------------
// Instantiations.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class GradyKippScalarDamage<Dim<1> >;
  template class GradyKippScalarDamage<Dim<2> >;
  template class GradyKippScalarDamage<Dim<3> >;
}
