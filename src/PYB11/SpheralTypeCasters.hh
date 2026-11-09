//------------------------------------------------------------------------------
// https://pybind11.readthedocs.io/en/stable/advanced/cast/stl.html#making-opaque-types
// The documentation requires opaque declarations before use in every relevant
// translation unit, typically through a common header.
//
// https://github.com/pybind/pybind11/issues/1055
// The discussion explains ODR violations caused by inconsistent type_caster
// definitions across translation units.
//
// https://github.com/pybind/pybind11/issues/1940
// https://stackoverflow.com/questions/58169847/mixing-type-conversions-and-opaque-types-with-pybind11
// The discussions show implicitly_convertible<py::list, std::vector<T>>()
// restoring Python-list inputs for opaque, bound vectors.
//
// https://pybind11.readthedocs.io/en/stable/advanced/classes.html#implicit-conversions
// The documentation describes implicit-conversion registration and requires
// the destination vector class to be bound before registration.
//
// https://pybind11.readthedocs.io/en/stable/advanced/cast/custom.html
// The documentation describes type_caster specializations and load().
// The PythonVectorInput specialization uses load() to filter input categories;
// the bound vector's constructor performs the element conversion.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Shared opaque-container policy for all Spheral Python binding translation units.
//
// PYBIND11_MAKE_OPAQUE defines the corresponding pybind11::detail::type_caster
// specialization. Do not repeat those definitions in a generated module header.
// The generator's PYB11_SHARED_TYPE_CASTER_HEADER hook selects this file instead.
//
// This is the union of the current PYB11opaque lists and opaque=True STL bindings.
// Add new opaque container types here when adding bindings. All dimensions are
// declared because Geometry and the common PYB11opaque lists include all three,
// even when only a subset of dimensional binding classes is enabled.
//
// The shared policy defines one bound-object type_caster specialization per type.
// After each vector class is bound, the generator calls the helper below to accept
// the input kinds supported by pybind11's STL vector type_caster, using the class's
// iterable constructor. Results retain bound-object conversion behavior;
// list snapshots require separate API adapters.
//------------------------------------------------------------------------------
#ifndef SPHERAL_PYB11_TYPE_CASTERS_HH
#define SPHERAL_PYB11_TYPE_CASTERS_HH

// Python must precede system headers that can define _POSIX_C_SOURCE.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <map>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Geometry/Dimension.hh"
#include "Geometry/PolyClipperUtilities.hh"

// Container type_caster declarations do not need the full element definitions.
// Avoid including headers with Python conversion code before this policy.
namespace Spheral {
template<typename Dimension> class GeomPlane;
template<typename Dimension> struct RKCoefficients;
struct CellFaceFlag;
template<typename Dimension> class NodeList;
template<typename Dimension> class FluidNodeList;
template<typename Dimension> class SolidNodeList;
template<typename Dimension> class Physics;
template<typename Dimension> class Boundary;
template<typename Dimension> struct DomainNode;
template<typename Dimension> class GridCellIndex;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
}

// CXXTypes, SpheralCommon, Silo, polytope, and Solvers.
PYBIND11_MAKE_OPAQUE(std::vector<char>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<char>>);
PYBIND11_MAKE_OPAQUE(std::vector<unsigned>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<unsigned>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::uint64_t>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::uint64_t>>);
PYBIND11_MAKE_OPAQUE(std::vector<int>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<int>>);
PYBIND11_MAKE_OPAQUE(std::vector<float>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<float>>);
PYBIND11_MAKE_OPAQUE(std::vector<double>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<double>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::string>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<std::string>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::pair<double, double>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::pair<double, std::string>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::pair<unsigned, unsigned>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::pair<std::uint64_t, std::uint64_t>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::pair<std::string, std::string>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::tuple<int, int, int>>);
PYBIND11_MAKE_OPAQUE(std::map<std::string, double>);
PYBIND11_MAKE_OPAQUE(std::map<int, std::string>);

// Geometry and RK value containers, including both PolyClipper adapters.
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::Vector>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::Vector>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::Tensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::Tensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::SymTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::SymTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::ThirdRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::ThirdRankTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::FourthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::FourthRankTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::FifthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::FifthRankTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::FacetedVolume>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::FacetedVolume>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomPlane<Spheral::Dim<1>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::RKCoefficients<Spheral::Dim<1>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::Vector>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::Vector>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::Tensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::Tensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::SymTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::SymTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::ThirdRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::ThirdRankTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::FourthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::FourthRankTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::FifthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::FifthRankTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::FacetedVolume>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::FacetedVolume>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomPlane<Spheral::Dim<2>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::RKCoefficients<Spheral::Dim<2>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::Vector>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::Vector>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::Tensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::Tensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::SymTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::SymTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::ThirdRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::ThirdRankTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::FourthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::FourthRankTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::FifthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::FifthRankTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::FacetedVolume>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::FacetedVolume>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomPlane<Spheral::Dim<3>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::RKCoefficients<Spheral::Dim<3>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomFacet2d>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomFacet3d>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::CellFaceFlag>);
PYBIND11_MAKE_OPAQUE(std::vector<PolyClipper::Vertex2d<>>);
PYBIND11_MAKE_OPAQUE(std::vector<PolyClipper::Vertex3d<>>);
PYBIND11_MAKE_OPAQUE(std::vector<PolyClipper::Vertex2d<Spheral::GeomVectorAdapter<2>>>);
PYBIND11_MAKE_OPAQUE(std::vector<PolyClipper::Vertex3d<Spheral::GeomVectorAdapter<3>>>);

// NodeList, Physics, Boundary, Utilities, and Neighbor containers.
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::NodeList<Spheral::Dim<1>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FluidNodeList<Spheral::Dim<1>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::SolidNodeList<Spheral::Dim<1>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Physics<Spheral::Dim<1>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Boundary<Spheral::Dim<1>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<std::pair<Spheral::NodeList<Spheral::Dim<1>>*, std::string>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::DomainNode<Spheral::Dim<1>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GridCellIndex<Spheral::Dim<1>>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::GridCellIndex<Spheral::Dim<1>>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::NodeList<Spheral::Dim<2>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FluidNodeList<Spheral::Dim<2>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::SolidNodeList<Spheral::Dim<2>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Physics<Spheral::Dim<2>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Boundary<Spheral::Dim<2>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<std::pair<Spheral::NodeList<Spheral::Dim<2>>*, std::string>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::DomainNode<Spheral::Dim<2>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GridCellIndex<Spheral::Dim<2>>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::GridCellIndex<Spheral::Dim<2>>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::NodeList<Spheral::Dim<3>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FluidNodeList<Spheral::Dim<3>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::SolidNodeList<Spheral::Dim<3>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Physics<Spheral::Dim<3>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Boundary<Spheral::Dim<3>>*>);
PYBIND11_MAKE_OPAQUE(std::vector<std::pair<Spheral::NodeList<Spheral::Dim<3>>*, std::string>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::DomainNode<Spheral::Dim<3>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GridCellIndex<Spheral::Dim<3>>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::GridCellIndex<Spheral::Dim<3>>>>);

// Field/FieldList value, pointer, and nested containers (also DataBase/FileIO).
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<1>, int>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<1>, int>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<1>, int>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<1>, double>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<1>, double>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<1>, double>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::Vector>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::Vector>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::Vector>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<1>, int>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<1>, int>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<1>, int>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<1>, double>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<1>, double>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<1>, double>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Vector>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::Tensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<1>, Spheral::Dim<1>::SymTensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<2>, int>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<2>, int>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<2>, int>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<2>, double>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<2>, double>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<2>, double>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::Vector>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::Vector>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::Vector>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<2>, int>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<2>, int>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<2>, int>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<2>, double>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<2>, double>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<2>, double>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Vector>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::Tensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<2>, Spheral::Dim<2>::SymTensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<3>, int>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<3>, int>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<3>, int>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<3>, double>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<3>, double>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<3>, double>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::Vector>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::Vector>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::Vector>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Field<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<3>, int>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<3>, int>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<3>, int>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<3>, double>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<3>, double>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<3>, double>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Vector>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::Tensor>*>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::FieldList<Spheral::Dim<3>, Spheral::Dim<3>::SymTensor>*>);

namespace Spheral {

// A source-filter marker for implicitly_convertible, not a Python-visible class
// or a replacement for any native vector type. It holds no converted data.
struct PythonVectorInput {};

}

PYBIND11_NAMESPACE_BEGIN(PYBIND11_NAMESPACE)
namespace detail {

template<>
struct type_caster<Spheral::PythonVectorInput> {
  PYBIND11_TYPE_CASTER(Spheral::PythonVectorInput, const_name("STL vector input"));

  bool load(handle source, bool /* convert */) {
    // Use the same input-category check as stl.h's list_caster. In particular,
    // bound vectors qualify as sequences, but strings/bytes/plain dicts do not.
    // This is a pybind11 internal helper: review it when updating pybind11.
    return source && object_is_convertible_to_std_vector(source);
  }
};

} // namespace detail
PYBIND11_NAMESPACE_END(PYBIND11_NAMESPACE)

namespace Spheral {

// Run during module initialization, once after each vector class is registered.
// Merely including this header must not register conversions: the Python class
// does not exist yet. The source filter checks the kind of Python input; the
// existing vector constructor checks its elements. Matching bound-vector inputs
// keep their direct conversion path; other accepted inputs create temporaries.
template<typename Vector>
inline void registerPythonVectorConversion() {
  pybind11::implicitly_convertible<PythonVectorInput, Vector>();
}

}

// The generator invokes this template immediately after py::bind_vector.
// An object-like macro keeps vector element types containing commas intact.
#define PYB11_VECTOR_BINDING_HOOK ::Spheral::registerPythonVectorConversion

#endif
