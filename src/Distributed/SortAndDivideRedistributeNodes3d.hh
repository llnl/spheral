//---------------------------------Spheral++----------------------------------//
// SortAndDivideRedistributeNodes3d -- 3-D implementation of the sort and 
// divide algorithm for domain decomposition.
//
// Created by JMO, Sun Dec  5 21:44:28 PST 2004
//----------------------------------------------------------------------------//
#ifndef Spheral_SortAndDivideRedistributeNodes3d_hh
#define Spheral_SortAndDivideRedistributeNodes3d_hh

#include "SortAndDivideRedistributeNodes.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class NodeList;
template<typename Dimension> class Boundary;
template<typename Dimension, typename DataType> class FieldList;

class SortAndDivideRedistributeNodes3d: public SortAndDivideRedistributeNodes<Dim<3> > {

public:
  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<3>;
  using Scalar = Dim<3>::Scalar;
  using Vector = Dim<3>::Vector;
  using Tensor = Dim<3>::Tensor;
  using SymTensor = Dim<3>::SymTensor;

  // Constructors
  SortAndDivideRedistributeNodes3d(const double HExtent);

  // Destructor
  virtual ~SortAndDivideRedistributeNodes3d() = default;

  // Given a Spheral++ data base of NodeLists, repartition it among the processors.
  virtual void redistributeNodes(DataBase<Dim<3> >& dataBase,
                                 std::vector<Boundary<Dim<3> >*> boundaries = std::vector<Boundary<Dim<3> >*>());

  // Compute the appropriate number of domains in each dimension for a given
  // shape tensor EigenStruct.
  std::vector< std::vector<int> > domainsPerChunk(const SymTensor::EigenStructType& shapeTensor) const;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copy, or assignment operations.
  SortAndDivideRedistributeNodes3d();
  SortAndDivideRedistributeNodes3d(const SortAndDivideRedistributeNodes3d&);
  SortAndDivideRedistributeNodes3d& operator=(const SortAndDivideRedistributeNodes3d&);

};

}

#endif
