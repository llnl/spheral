text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Utilities/nodeBoundingBoxes.cc"

namespace Spheral {

  template Field<Dim<%(ndim)s>, std::pair<Dim<%(ndim)s>::Vector, Dim<%(ndim)s>::Vector>> nodeBoundingBoxes(const NodeList<Dim<%(ndim)s>>&, 
                                                                                                           const Dim<%(ndim)s>::Scalar);
  template FieldList<Dim<%(ndim)s>, std::pair<Dim<%(ndim)s>::Vector, Dim<%(ndim)s>::Vector>> nodeBoundingBoxes(const DataBase<Dim<%(ndim)s>>& dataBase);

}
"""
