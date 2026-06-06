//---------------------------------Spheral++----------------------------------//
// RigidBoundary -- Apply a Reflecting boundary condition to Spheral++
// Fields.
//
// Created by JMO, Wed Feb 16 21:01:02 PST 2000
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "Geometry/GeomPlane.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/DBC.hh"

#include "RigidBoundary.hh"

using std::vector;

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
RigidBoundary<Dimension>::RigidBoundary():
  ReflectingBoundary<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given plane.
//------------------------------------------------------------------------------
template<typename Dimension>
RigidBoundary<Dimension>::
RigidBoundary(const GeomPlane<Dimension>& plane):
  ReflectingBoundary<Dimension>(plane) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
RigidBoundary<Dimension>::~RigidBoundary() {
}

//------------------------------------------------------------------------------
// Apply the ghost boundary condition to fields of different DataTypes.
//------------------------------------------------------------------------------
// Specialization for Vector fields.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  REQUIRE(this->valid());
  if (field.name() == HydroFieldNames::velocity) {
    ReflectingBoundary<Dimension>::applyGhostBoundary(field);
  } else {
    Boundary<Dimension>::applyGhostBoundary(dynamic_cast<FieldBase<Dimension>&>(field));
  }
}

// Specialization for Tensor fields.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::Tensor>& field) const {
  Boundary<Dimension>::applyGhostBoundary(dynamic_cast<FieldBase<Dimension>&>(field));
}

// Specialization for symmetric tensors.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::SymTensor>& field) const {
  Boundary<Dimension>::applyGhostBoundary(dynamic_cast<FieldBase<Dimension>&>(field));
}

// Specialization for third rank tensors.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>& field) const {
  Boundary<Dimension>::applyGhostBoundary(dynamic_cast<FieldBase<Dimension>&>(field));
}

// Specialization for fourth rank tensors.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FourthRankTensor>& field) const {
  Boundary<Dimension>::applyGhostBoundary(dynamic_cast<FieldBase<Dimension>&>(field));
}

// Specialization for fifth rank tensors.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FifthRankTensor>& field) const {
  Boundary<Dimension>::applyGhostBoundary(dynamic_cast<FieldBase<Dimension>&>(field));
}

// Specialization for FacetedVolume.
template<typename Dimension>
void
RigidBoundary<Dimension>::
applyGhostBoundary(Field<Dimension, typename Dimension::FacetedVolume>& field) const {
  Boundary<Dimension>::applyGhostBoundary(dynamic_cast<FieldBase<Dimension>&>(field));
}

//------------------------------------------------------------------------------
// Enforce the boundary condition on the set of nodes in violation of the 
// boundary.
//------------------------------------------------------------------------------
// Specialization for vector fields.  Apply the reflection operator to x and v.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Vector>& field) const {
  REQUIRE(this->valid());
  if (field.name() == HydroFieldNames::velocity) ReflectingBoundary<Dimension>::applyGhostBoundary(field);
}

// Specialization for tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::Tensor>&) const {
}

// Specialization for tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::SymTensor>&) const {
}

// Specialization for third rank tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::ThirdRankTensor>&) const {
}

// Specialization for fourth rank tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FourthRankTensor>&) const {
}

// Specialization for fifth rank tensor fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FifthRankTensor>&) const {
}

// Specialization for FacetedVolume fields.  No-op.
template<typename Dimension>
void
RigidBoundary<Dimension>::
enforceBoundary(Field<Dimension, typename Dimension::FacetedVolume>&) const {
}

}
