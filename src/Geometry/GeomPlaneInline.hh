namespace Spheral {

//------------------------------------------------------------------------------
// Global comparator functions.
//------------------------------------------------------------------------------
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator==(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane == vec;
}
  
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator!=(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane != vec;
}
  
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator<(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane > vec;
}
  
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator>(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane < vec;
}
  
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator<=(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane >= vec;
}
  
template<typename Dimension>
SPHERAL_HOST_DEVICE 
inline
bool
operator>=(const typename Dimension::Vector& vec, const GeomPlane<Dimension>& plane) {
  return plane <= vec;
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::istream&
operator>>(std::istream& is, GeomPlane<Dimension>& plane) {
  std::string parenthesis;
  typename Dimension::Vector point;
  typename Dimension::Vector normal;
  is >> parenthesis >> point >> normal >> parenthesis;
  plane.point(point);
  plane.normal(normal);
  return is;
}

//------------------------------------------------------------------------------
// Output (std::ostream) operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::ostream&
operator<<(std::ostream& os, const GeomPlane<Dimension>& plane) {
  os << "( " << plane.point() << plane.normal() << ")";
  return os;
}

}
