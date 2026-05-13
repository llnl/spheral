//---------------------------------Spheral++----------------------------------//
// GeometryScaling
//
// Free functions for scaling/unscaling FieldLists by the geometric factor
// appropriate to the coordinate system (2πr for RZ, r² for Spherical).
// These replace the hand-written loops scattered through RZ variant classes.
//
// The ghost parameter controls whether ghost nodes are included in the
// scaling loop (true = all nodes, false = internal only).
//----------------------------------------------------------------------------//
#ifndef __Spheral_GeometryScaling__
#define __Spheral_GeometryScaling__

#include "Field/FieldList.hh"
#include "Geometry/Dimension.hh"
#include "Geometry/GeometryRegistrar.hh"
#include "Utilities/DBC.hh"

#include <cmath>
#include <tuple>

namespace Spheral {

// So that the scaling is exactly reversible
inline double safeRadius(double r) {
  constexpr double fuzz = 1.0e-30;
  return std::max(fuzz, std::abs(r));
}

//------------------------------------------------------------------------------
// Dim<1>: scale by r² for spherical coordinates.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
scaleForGeometry(const FieldList<Dim<1>, Dim<1>::Vector>& position,
                 FieldList<Dim<1>, DataType>& field,
                 const bool ghost = true) {
  if (GeometryRegistrar::coords() == CoordinateType::Spherical) {
    for (auto nli = 0u; nli < field.numFields(); ++nli) {
      const auto n = ghost ? field[nli]->numElements() : field[nli]->numInternalElements();
      for (auto ni = 0u; ni < n; ++ni) {
        const auto r = safeRadius(position(nli, ni).x());
        field(nli, ni) *= r * r;
      }
    }
  }
}

template<typename DataType>
inline
void
unscaleFromGeometry(const FieldList<Dim<1>, Dim<1>::Vector>& position,
                    FieldList<Dim<1>, DataType>& field,
                    const bool ghost = true) {
  if (GeometryRegistrar::coords() == CoordinateType::Spherical) {
    for (auto nli = 0u; nli < field.numFields(); ++nli) {
      const auto n = ghost ? field[nli]->numElements() : field[nli]->numInternalElements();
      for (auto ni = 0u; ni < n; ++ni) {
        const auto r = safeRadius(position(nli, ni).x());
        field(nli, ni) /= r * r;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Dim<2>: scale by 2πr for RZ (cylindrical) coordinates.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
scaleForGeometry(const FieldList<Dim<2>, Dim<2>::Vector>& position,
                 FieldList<Dim<2>, DataType>& field,
                 const bool ghost = true) {
  if (GeometryRegistrar::coords() == CoordinateType::RZ) {
    for (auto nli = 0u; nli < field.numFields(); ++nli) {
      const auto n = ghost ? field[nli]->numElements() : field[nli]->numInternalElements();
      for (auto ni = 0u; ni < n; ++ni) {
        const auto r = safeRadius(position(nli, ni).y());
        field(nli, ni) *= 2.0 * M_PI * r;
      }
    }
  }
}

template<typename DataType>
inline
void
unscaleFromGeometry(const FieldList<Dim<2>, Dim<2>::Vector>& position,
                    FieldList<Dim<2>, DataType>& field,
                    const bool ghost = true) {
  if (GeometryRegistrar::coords() == CoordinateType::RZ) {
    for (auto nli = 0u; nli < field.numFields(); ++nli) {
      const auto n = ghost ? field[nli]->numElements() : field[nli]->numInternalElements();
      for (auto ni = 0u; ni < n; ++ni) {
        const auto r = safeRadius(position(nli, ni).y());
        field(nli, ni) /= 2.0 * M_PI * r;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Dim<3>: no geometric scaling needed.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
scaleForGeometry(const FieldList<Dim<3>, Dim<3>::Vector>& /*position*/,
                 FieldList<Dim<3>, DataType>& /*field*/,
                 const bool /*ghost*/ = true) {
}

template<typename DataType>
inline
void
unscaleFromGeometry(const FieldList<Dim<3>, Dim<3>::Vector>& /*position*/,
                    FieldList<Dim<3>, DataType>& /*field*/,
                    const bool /*ghost*/ = true) {
}

//------------------------------------------------------------------------------
// 3-arg overloads: copy and then apply scaling
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
scaleForGeometry(const FieldList<Dimension, typename Dimension::Vector>& position,
                 const FieldList<Dimension, DataType>& field,
                 FieldList<Dimension, DataType>& output,
                 const bool ghost = true) {
  output.assignFields(field);
  scaleForGeometry(position, output, ghost);
}

template<typename Dimension, typename DataType>
inline
void
unscaleFromGeometry(const FieldList<Dimension, typename Dimension::Vector>& position,
                    const FieldList<Dimension, DataType>& field,
                    FieldList<Dimension, DataType>& output,
                    const bool ghost = true) {
  output.assignFields(field);
  unscaleFromGeometry(position, output, ghost);
}

//------------------------------------------------------------------------------
// RegionScaling — single RAII guard for geometry scaling regions.
//   Inverse=true:  unscale on construction, rescale on destruction.
//   Inverse=false: scale on construction, unscale on destruction.
//   Ghost=true:    operate on all nodes (internal + ghost).
//   Ghost=false:   operate on internal nodes only.
//------------------------------------------------------------------------------
template<bool Inverse, bool Ghost, typename Dimension, typename... Fields>
class [[nodiscard]] RegionScaling {
public:
  RegionScaling(const FieldList<Dimension, typename Dimension::Vector>& pos,
                Fields&... fields)
    : mPosition(pos), mFields(fields...) {
    if constexpr (Inverse)
      std::apply([&](auto&... fs) { (unscaleFromGeometry(mPosition, fs, Ghost), ...); }, mFields);
    else
      std::apply([&](auto&... fs) { (scaleForGeometry(mPosition, fs, Ghost), ...); }, mFields);
  }
  ~RegionScaling() {
    if constexpr (Inverse)
      std::apply([&](auto&... fs) { (scaleForGeometry(mPosition, fs, Ghost), ...); }, mFields);
    else
      std::apply([&](auto&... fs) { (unscaleFromGeometry(mPosition, fs, Ghost), ...); }, mFields);
  }
  RegionScaling(const RegionScaling&) = delete;
  RegionScaling& operator=(const RegionScaling&) = delete;
  bool applied() { return GeometryRegistrar::coords() != CoordinateType::Cartesian; }

private:
  const FieldList<Dimension, typename Dimension::Vector>& mPosition;
  std::tuple<Fields&...> mFields;
};

// Factory functions — names match the old classes, lowercase to signal functions.
// C++17 guaranteed copy elision means no copy/move needed.
template<typename Dimension, typename... Fields>
auto unscaledRegion(const FieldList<Dimension, typename Dimension::Vector>& pos, Fields&... fields) {
  return RegionScaling<true, true, Dimension, Fields...>(pos, fields...);
}

template<typename Dimension, typename... Fields>
auto unscaledInternalRegion(const FieldList<Dimension, typename Dimension::Vector>& pos, Fields&... fields) {
  return RegionScaling<true, false, Dimension, Fields...>(pos, fields...);
}

template<typename Dimension, typename... Fields>
auto scaledRegion(const FieldList<Dimension, typename Dimension::Vector>& pos, Fields&... fields) {
  return RegionScaling<false, true, Dimension, Fields...>(pos, fields...);
}

template<typename Dimension, typename... Fields>
auto scaledInternalRegion(const FieldList<Dimension, typename Dimension::Vector>& pos, Fields&... fields) {
  return RegionScaling<false, false, Dimension, Fields...>(pos, fields...);
}

} // end namespace Spheral

#endif
