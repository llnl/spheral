text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SVPH/computeSVPHCorrectionsOnFaces.cc"
#include "Geometry/Dimension.hh"
#include "Physics/Physics.hh"
#include "Boundary/Boundary.hh"

namespace Spheral {

template 
void
computeSVPHCorrectionsOnFaces(const Mesh<Dim< %(ndim)s > >&,
                              const DataBase<Dim< %(ndim)s > >&,
                              const TableKernel<Dim< %(ndim)s > >&,
                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>&,
                              const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::SymTensor>&,
                              const Physics<Dim< %(ndim)s > >::ConstBoundaryIterator&,
                              const Physics<Dim< %(ndim)s > >::ConstBoundaryIterator&,
                              std::vector<Dim< %(ndim)s >::Scalar>&,
                              std::vector<Dim< %(ndim)s >::Vector>&);

}
"""
