//---------------------------------Spheral++----------------------------------//
// FieldListSet -- A container class for FieldLists.  Provides a convenient
//   container for a set of associated FieldLists, allowing us to pass them
//   in bundles through routines such as FieldListFunctions.
//
// Created by JMO, Wed Jun  5 21:47:20 PDT 2002
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldListSet_hh__
#define __Spheral_FieldListSet_hh__

#include "FieldList.hh"
#include <vector>

namespace Spheral {

template<typename Dimension>
class FieldListSet {
public:
  std::vector< FieldList<Dimension, typename Dimension::Scalar> > ScalarFieldLists;
  std::vector< FieldList<Dimension, typename Dimension::Vector> > VectorFieldLists;
  std::vector< FieldList<Dimension, typename Dimension::Tensor> > TensorFieldLists;
  std::vector< FieldList<Dimension, typename Dimension::SymTensor> > SymTensorFieldLists;
};

}

#endif
