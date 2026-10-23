//---------------------------------Spheral++----------------------------------//
// Neighbor -- An abstract base class for the Neighbor objects.
//
// Created by J. Michael Owen, Sun Nov 12 10:33:55 2000
// Split into Value/View pattern Tue Sep 22 13:55:08 PDT 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_Neighbor__
#define __Spheral_Neighbor__

#include "Neighbor/NeighborSearchType.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <utility>     // std::pair
#include <span>

namespace Spheral {

// template<typename Dimension> class DataBase;
// template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension>
class Neighbor {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using NodeID = std::pair<size_t, size_t>;  // (NodeList ID, node ID)

  // Constructors and destructors
  Neighbor(DataBase<Dimension>& dataBase,
           const NeighborSearchType searchType,
           const Scalar kernelExtent);
  Neighbor(const Neighbor& rhs) = default;
  Neighbor(Neighbor&& rhs) = default;
  virtual ~Neighbor() = default;

  // Public access to our basic state data
  const FieldList<Dimension, Vector>& nodeExtent()       const { return mNodeExtent; }

  // Set internal data
  void neighborSearchType(const NeighborSearchType& x)         { mSearchType = x; }
  void kernelExtent(const Scalar x)                            { VERIFY(x > 0.0); mKernelExtent = x; }

  //.............................................................................
  // Requirements for rebuilding the internal data of a Neighbor object.
  // Expected to be CPU Value object tasks for now.
  // Force the update of internal data for the NodeList.
  virtual void updateNodes() = 0;
  virtual void updateNodes(const std::span<NodeID>& nodeIDs) = 0;

  // Reinitialize based on the desired target h and box dimensions.
  virtual void reinitialize() = 0;
  virtual void reinitialize(const Vector& xmin,
                            const Vector& xmax,
                            const Scalar htarget) = 0;

protected:
  //-------------------------- Protected Interface --------------------------//
  DataBase<Dimension>& mDataBase;
  FieldList<Dimension, Vector> mNodeExtent;
  NeighborSearchType mSearchType;
  Scalar mKernelExtent;
};

}

#endif
