//---------------------------------Spheral++----------------------------------//
// PeriodicPlanarBoundary -- The nested class member of Periodic Boundary that
// does the actual work.
//
// Created by JMO, Wed Apr 19 15:00:50 PDT 2000
//----------------------------------------------------------------------------//

// #include "PeriodicBoundary.hh"
#include "Utilities/DBC.hh"

namespace {

//------------------------------------------------------------------------------
// Map a faceted volume
//------------------------------------------------------------------------------
template<typename BC, typename Dimension>
struct MapFacetedVolume {
  typename Dimension::FacetedVolume 
  static doit(const BC& bc,
       const typename Dimension::FacetedVolume& poly) {
    typedef typename Dimension::Vector Vector;
    typedef typename Dimension::FacetedVolume FacetedVolume;
    const auto& enterPlane = bc.enterPlane();
    const auto& exitPlane = bc.exitPlane();
    const auto& verts0 = poly.vertices();
    const auto& facets = poly.facetVertices();
    vector<Vector> verts1;
    for (auto& vert: verts0) verts1.push_back(bc.mapPosition(vert, enterPlane, exitPlane));
    return FacetedVolume(verts1, facets);
  }
};

template<typename BC>
struct MapFacetedVolume<BC, Dim<1>> {
  Dim<1>::FacetedVolume 
  static doit(const BC& bc,
       const Dim<1>::FacetedVolume& poly) {
    const auto& enterPlane = bc.enterPlane();
    const auto& exitPlane = bc.exitPlane();
    return Dim<1>::FacetedVolume(bc.mapPosition(poly.center(),
                                                enterPlane,
                                                exitPlane),
                                 poly.extent());
  }
};

}

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
PeriodicPlanarBoundary():
  PlanarBoundary<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given planes.
//------------------------------------------------------------------------------
template<typename Dimension>
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
PeriodicPlanarBoundary(const GeomPlane<Dimension>& plane1,
                       const GeomPlane<Dimension>& plane2):
  PlanarBoundary<Dimension>(plane1, plane2) {

  // Once we're done the boundary condition should be in a valid state.
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
~PeriodicPlanarBoundary() {
}

//------------------------------------------------------------------------------
// Apply the boundary condition to the ghost values in the given field.
// FieldBase forward
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(FieldBase<Dimension>& fieldBase) const {
  Boundary<Dimension>::applyGhostBoundary(fieldBase);
}

// FacetedVolume fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
applyGhostBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {

  REQUIRE(valid());

  // Apply the boundary condition to all the ghost node values.
  const NodeList<Dimension>& nodeList = field.nodeList();
  CHECK(this->controlNodes(nodeList).size() == this->ghostNodes(nodeList).size());
  auto controlItr = this->controlBegin(nodeList);
  auto ghostItr = this->ghostBegin(nodeList);
  for (; controlItr < this->controlEnd(nodeList); ++controlItr, ++ghostItr) {
    CHECK(ghostItr < this->ghostEnd(nodeList));
    CHECK(*controlItr < nodeList.numNodes());
    CHECK(*ghostItr >= nodeList.firstGhostNode() and *ghostItr < nodeList.numNodes());
    field(*ghostItr) = MapFacetedVolume<PeriodicBoundary<Dimension>::PeriodicPlanarBoundary, Dimension>::doit(*this, field(*controlItr));
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary condition on nodes in violation for the given field.
//------------------------------------------------------------------------------
// FacetedVolume fields.
template<typename Dimension>
void
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::
enforceBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  REQUIRE(valid());
  const auto& nodeList = field.nodeList();
  const auto& violationNodes = this->violationNodes(nodeList);
  for (auto i: violationNodes) {
    CHECK(i < nodeList.numInternalNodes());
    field(i) = MapFacetedVolume<PeriodicBoundary<Dimension>::PeriodicPlanarBoundary, Dimension>::doit(*this, field(i));
  }
}

//------------------------------------------------------------------------------
// Test if the periodic planar boundary is minimally valid.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PeriodicBoundary<Dimension>::PeriodicPlanarBoundary::valid() const {
  return (PlanarBoundary<Dimension>::valid());
}
