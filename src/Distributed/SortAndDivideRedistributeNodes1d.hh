//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes1d -- 1-D implementation of the sort and 
// divide algorithm for domain decomposition.
//
// Created by JMO, Thu Dec 2 11:01:53 2004
//----------------------------------------------------------------------------//
#ifndef Spheral_SortAndDivideRedistributeNodes1d_hh
#define Spheral_SortAndDivideRedistributeNodes1d_hh

#include "SortAndDivideRedistributeNodes.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class NodeList;
template<typename Dimension> class Boundary;
template<typename Dimension, typename DataType> class FieldList;

class SortAndDivideRedistributeNodes1d: public SortAndDivideRedistributeNodes<Dim<1> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<1>;
  using Scalar = Dim<1>::Scalar;
  using Vector = Dim<1>::Vector;
  using Tensor = Dim<1>::Tensor;
  using SymTensor = Dim<1>::SymTensor;

  // Constructors
  SortAndDivideRedistributeNodes1d() = default;
  SortAndDivideRedistributeNodes1d(const SortAndDivideRedistributeNodes1d&) = delete;
  SortAndDivideRedistributeNodes1d& operator=(const SortAndDivideRedistributeNodes1d&) = delete;
  virtual ~SortAndDivideRedistributeNodes1d() = default;

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  virtual void redistributeNodes(DataBase<Dim<1> >& dataBase,
                                 std::vector<Boundary<Dim<1> >*> boundaries = std::vector<Boundary<Dim<1> >*>());
};

}

#endif
