#include <boost/geometry.hpp>

#include "Geometry/Dimension.hh"

//------------------------------------------------------------------------------
// GeomVector<nDim> -> Boost.Geometry
//------------------------------------------------------------------------------
namespace boost {
namespace geometry {
namespace traits {

// Adapt Spheral::GeomVector<2> to Boost.Geometry

template<int nDim> struct tag<Spheral::GeomVector<nDim>>                                { using type = point_tag; };
template<int nDim> struct dimension<Spheral::GeomVector<nDim>> : boost::mpl::int_<nDim> {};
template<int nDim> struct coordinate_type<Spheral::GeomVector<nDim>>                    { using type = double; };
template<int nDim> struct coordinate_system<Spheral::GeomVector<nDim>>                  { using type = cs::cartesian; };

template<int nDim, std::size_t Index>
struct access<Spheral::GeomVector<nDim>, Index> {
  static_assert(Index < nDim, "Index out of dimensional range");
  using Point = Spheral::GeomVector<nDim>;
  using CoordinateType = typename coordinate_type<Point>::type;
  static inline CoordinateType get(Point const& p) { return p[Index]; }
  static inline void set(Point& p, CoordinateType const& value) { p[Index] = value; }
};

}
}
} // namespace boost::geometry::traits

