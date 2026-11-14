//------------------------------------------------------------------------------
// Modified ideas from CRKSPH using the SVPH normalized kernel estimates.
// This version computes the corrections on the faces of a mesh, rather than on
// the points themselves.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeSVPHCorrectionsOnFaces__
#define __Spheral__computeSVPHCorrectionsOnFaces__

#include "Geometry/Dimension.hh"

#include <vector>

namespace Spheral {

// Forward declarations.
template<typename Dimension> class Mesh;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;

template<typename Dimension, typename BoundaryIterator>
void
computeSVPHCorrectionsOnFaces(const Mesh<Dimension>& mesh,
                              const DataBase<Dimension>& dataBase,
                              const TableKernel<Dimension>& W,
                              const FieldList<Dimension, typename Dimension::Scalar>& volume,
                              const FieldList<Dimension, typename Dimension::Vector>& position,
                              const FieldList<Dimension, typename Dimension::SymTensor>& H,
                              const BoundaryIterator& boundaryBegin,
                              const BoundaryIterator& boundaryEnd,
                              std::vector<typename Dimension::Scalar>& A,
                              std::vector<typename Dimension::Vector>& B);

}

#endif
