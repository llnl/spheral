#include "Geometry/Dimension.hh"
#include "GridCellPlane.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Define a global function which returns a vector of GridCellIndex objects
// within a given range of grid cells.
//------------------------------------------------------------------------------
std::vector<GridCellIndex<Dim<1> > >
GridCellIndexRange(const GridCellIndex<Dim<1> >& gridCellMin,
                   const GridCellIndex<Dim<1> >& gridCellMax);
std::vector<GridCellIndex<Dim<2> > >
GridCellIndexRange(const GridCellIndex<Dim<2> >& gridCellMin,
                   const GridCellIndex<Dim<2> >& gridCellMax);
std::vector<GridCellIndex<Dim<3> > >
GridCellIndexRange(const GridCellIndex<Dim<3> >& gridCellMin,
                   const GridCellIndex<Dim<3> >& gridCellMax);

//------------------------------------------------------------------------------
// Similarly, define an iterator-like function that will sequentially return
// the next grid cell in a sequence, given the current grid cell, the min, and 
// the max.
//------------------------------------------------------------------------------
inline
void
incrementGridCell(GridCellIndex<Dim<1> >& gridCell,
                  const GridCellIndex<Dim<1> >& gridCellMin,
                  const GridCellIndex<Dim<1> >& gridCellMax) {
  CONTRACT_VAR(gridCellMin);
  CONTRACT_VAR(gridCellMax);
  REQUIRE(gridCellMax.xIndex() >= gridCellMin.xIndex());
  REQUIRE(gridCell >= gridCellMin && gridCell <= gridCellMax);
  gridCell.xIndex(gridCell.xIndex() + 1);
}

inline
void
incrementGridCell(GridCellIndex<Dim<2> >& gridCell,
                  const GridCellIndex<Dim<2> >& gridCellMin,
                  const GridCellIndex<Dim<2> >& gridCellMax) {
  REQUIRE(gridCellMax.xIndex() >= gridCellMin.xIndex());
  REQUIRE(gridCellMax.yIndex() >= gridCellMin.yIndex());
  REQUIRE(gridCell >= gridCellMin && gridCell <= gridCellMax);
  gridCell.xIndex(gridCell.xIndex() + 1);
  if (gridCell.xIndex() > gridCellMax.xIndex()) {
    gridCell.xIndex(gridCellMin.xIndex());
    gridCell.yIndex(gridCell.yIndex() + 1);
  }
}

inline
void
incrementGridCell(GridCellIndex<Dim<3> >& gridCell,
                  const GridCellIndex<Dim<3> >& gridCellMin,
                  const GridCellIndex<Dim<3> >& gridCellMax) {
  REQUIRE(gridCellMax.xIndex() >= gridCellMin.xIndex());
  REQUIRE(gridCellMax.yIndex() >= gridCellMin.yIndex());
  REQUIRE(gridCellMax.zIndex() >= gridCellMin.zIndex());
  REQUIRE(gridCell >= gridCellMin && gridCell <= gridCellMax);
  gridCell.xIndex(gridCell.xIndex() + 1);
  if (gridCell.xIndex() > gridCellMax.xIndex()) {
    gridCell.xIndex(gridCellMin.xIndex());
    gridCell.yIndex(gridCell.yIndex() + 1);
    if (gridCell.yIndex() > gridCellMax.yIndex()) {
      gridCell.yIndex(gridCellMin.yIndex());
      gridCell.zIndex(gridCell.zIndex() + 1);
    }
  }
}

//------------------------------------------------------------------------------
// Multiply a GridCellIndex by a double, and return a Vector position.
//------------------------------------------------------------------------------
inline
Dim<1>::Vector
operator*(const GridCellIndex<Dim<1> >& lhs, const double rhs) {
  return Dim<1>::Vector(lhs.xIndex()*rhs);
}

inline
Dim<2>::Vector
operator*(const GridCellIndex<Dim<2> >& lhs, const double rhs) {
  return Dim<2>::Vector(lhs.xIndex()*rhs,
                        lhs.yIndex()*rhs);
}

inline
Dim<3>::Vector
operator*(const GridCellIndex<Dim<3> >& lhs, const double rhs) {
  return Dim<3>::Vector(lhs.xIndex()*rhs,
                        lhs.yIndex()*rhs,
                        lhs.zIndex()*rhs);
}

template<typename Dimension>
inline
typename Dimension::Vector
operator*(const double lhs, const GridCellIndex<Dimension>& rhs) {
  return rhs*lhs;
}

//------------------------------------------------------------------------------
// Empty Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GridCellIndex<Dimension>::GridCellIndex():
  GridCellIndexBase<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given indices.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>::GridCellIndex(int /*xIndex*/):
  GridCellIndexBase<Dimension>() {
  VERIFY(false);
}

template<typename Dimension>
inline
GridCellIndex<Dimension>::GridCellIndex(int /*xIndex*/, int /*yIndex*/):
  GridCellIndexBase<Dimension>() {
  VERIFY(false);
}

template<typename Dimension>
inline
GridCellIndex<Dimension>::GridCellIndex(int /*xIndex*/, int /*yIndex*/, int /*zIndex*/):
  GridCellIndexBase<Dimension>() {
  VERIFY(false);
}

template<>
inline
GridCellIndex<Dim<1> >::GridCellIndex(int xIndex):
  GridCellIndexBase<Dim<1> >(xIndex) {
}

template<>
inline
GridCellIndex<Dim<2> >::GridCellIndex(int xIndex, int yIndex):
  GridCellIndexBase<Dim<2> >(xIndex, yIndex) {
}

template<>
inline
GridCellIndex<Dim<3> >::GridCellIndex(int xIndex, int yIndex, int zIndex):
  GridCellIndexBase<Dim<3> >(xIndex, yIndex, zIndex) {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GridCellIndex<Dimension>::GridCellIndex(const GridCellIndex<Dimension>& rhs):
  GridCellIndexBase<Dimension>(rhs) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
GridCellIndex<Dimension>::~GridCellIndex() {
}

//------------------------------------------------------------------------------
// Return the X index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::xIndex() const {
  return this->mx;
}

//------------------------------------------------------------------------------
// Return the Y index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::yIndex() const {
  VERIFY(false);
  return 0;
}

template<>
inline
int
GridCellIndex<Dim<2> >::yIndex() const {
  return this->my;
}

template<>
inline
int
GridCellIndex<Dim<3> >::yIndex() const {
  return this->my;
}

//------------------------------------------------------------------------------
// Return the Z index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::zIndex() const {
  VERIFY(false);
  return 0;
}

template<>
inline
int
GridCellIndex<Dim<3> >::zIndex() const {
  return this->mz;
}

//------------------------------------------------------------------------------
// Set the X index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GridCellIndex<Dimension>::xIndex(int xIndex) {
  this->mx = xIndex;
}

//------------------------------------------------------------------------------
// Set the Y index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GridCellIndex<Dimension>::yIndex(int /*yIndex*/) {
  VERIFY(false);
}

template<>
inline
void
GridCellIndex<Dim<2> >::yIndex(int yIndex) {
  this->my = yIndex;
}

template<>
inline
void
GridCellIndex<Dim<3> >::yIndex(int yIndex) {
  this->my = yIndex;
}

//------------------------------------------------------------------------------
// Set the Z index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GridCellIndex<Dimension>::zIndex(int /*zIndex*/) {
  VERIFY(false);
}

template<>
inline
void
GridCellIndex<Dim<3> >::zIndex(int zIndex) {
  this->mz = zIndex;
}

//------------------------------------------------------------------------------
// Set the indices.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
GridCellIndex<Dimension>::setIndices(int /*xIndex*/) {
  VERIFY(false);
}

template<typename Dimension>
inline
void
GridCellIndex<Dimension>::setIndices(int /*xIndex*/,
                                     int /*yIndex*/) {
  VERIFY(false);
}

template<typename Dimension>
inline
void
GridCellIndex<Dimension>::setIndices(int /*xIndex*/,
                                     int /*yIndex*/,
                                     int /*zIndex*/) {
  VERIFY(false);
}

template<>
inline
void
GridCellIndex<Dim<1> >::setIndices(int xIndex) {
  this->mx = xIndex;
}

template<>
inline
void
GridCellIndex<Dim<2> >::setIndices(int xIndex, int yIndex) {
  this->mx = xIndex;
  this->my = yIndex;
}

template<>
inline
void
GridCellIndex<Dim<3> >::setIndices(int xIndex, int yIndex, int zIndex) {
  this->mx = xIndex;
  this->my = yIndex;
  this->mz = zIndex;
}

//------------------------------------------------------------------------------
// Access the elements by index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::operator()(int i) const {
  REQUIRE(i >= 0 and i < Dimension::nDim);
  return *(&(this->mx) + i);
}

template<typename Dimension>
inline
int&
GridCellIndex<Dimension>::operator()(int i) {
  REQUIRE(i >= 0 and i < Dimension::nDim);
  return *(&(this->mx) + i);
}

//------------------------------------------------------------------------------
// Assignment operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>&
GridCellIndex<Dimension>::operator=(const GridCellIndex<Dimension>& rhs) {
  GridCellIndexBase<Dimension>::operator=(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Provide the negative operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>
GridCellIndex<Dimension>::operator-() const {
  VERIFY(false);
  return *this;
}

template<>
inline
GridCellIndex<Dim<1> >
GridCellIndex<Dim<1> >::operator-() const {
  return GridCellIndex<Dim<1> >(-(this->mx));
}

template<>
inline
GridCellIndex<Dim<2> >
GridCellIndex<Dim<2> >::operator-() const {
  return GridCellIndex<Dim<2> >(-(this->mx),
                                -(this->my));
}

template<>
inline
GridCellIndex<Dim<3> >
GridCellIndex<Dim<3> >::operator-() const {
  return GridCellIndex<Dim<3> >(-(this->mx),
                                -(this->my),
                                -(this->mz));
}

//------------------------------------------------------------------------------
// Add a gridcell to another.
//------------------------------------------------------------------------------
template<>
inline
GridCellIndex<Dim<1> >
GridCellIndex<Dim<1> >::operator+(const GridCellIndex<Dim<1> >& rhs) const {
  return GridCellIndex<Dim<1> >(this->mx + rhs.mx);
}

template<>
inline
GridCellIndex<Dim<2> >
GridCellIndex<Dim<2> >::operator+(const GridCellIndex<Dim<2> >& rhs) const {
  return GridCellIndex<Dim<2> >(this->mx + rhs.mx,
                                this->my + rhs.my);
}

template<>
inline
GridCellIndex<Dim<3> >
GridCellIndex<Dim<3> >::operator+(const GridCellIndex<Dim<3> >& rhs) const {
  return GridCellIndex<Dim<3> >(this->mx + rhs.mx,
                                this->my + rhs.my,
                                this->mz + rhs.mz);
}

//------------------------------------------------------------------------------
// Subtract a gridcell from another.
//------------------------------------------------------------------------------
template<>
inline
GridCellIndex<Dim<1> >
GridCellIndex<Dim<1> >::operator-(const GridCellIndex<Dim<1> >& rhs) const {
  return GridCellIndex<Dim<1> >(this->mx - rhs.mx);
}

template<>
inline
GridCellIndex<Dim<2> >
GridCellIndex<Dim<2> >::operator-(const GridCellIndex<Dim<2> >& rhs) const {
  return GridCellIndex<Dim<2> >(this->mx - rhs.mx,
                                this->my - rhs.my);
}

template<>
inline
GridCellIndex<Dim<3> >
GridCellIndex<Dim<3> >::operator-(const GridCellIndex<Dim<3> >& rhs) const {
  return GridCellIndex<Dim<3> >(this->mx - rhs.mx,
                                this->my - rhs.my,
                                this->mz - rhs.mz);
}

//------------------------------------------------------------------------------
// Inline addition with a GridCellIndex.
//------------------------------------------------------------------------------
template<>
inline
GridCellIndex<Dim<1> >&
GridCellIndex<Dim<1> >::operator+=(const GridCellIndex<Dim<1> >& rhs) {
  this->mx += rhs.mx;
  return *this;
}

template<>
inline
GridCellIndex<Dim<2> >&
GridCellIndex<Dim<2> >::operator+=(const GridCellIndex<Dim<2> >& rhs) {
  this->mx += rhs.mx;
  this->my += rhs.my;
  return *this;
}

template<>
inline
GridCellIndex<Dim<3> >&
GridCellIndex<Dim<3> >::operator+=(const GridCellIndex<Dim<3> >& rhs) {
  this->mx += rhs.mx;
  this->my += rhs.my;
  this->mz += rhs.mz;
  return *this;
}

//------------------------------------------------------------------------------
// Inline subtraction with a GridCellIndex.
//------------------------------------------------------------------------------
template<>
inline
GridCellIndex<Dim<1> >&
GridCellIndex<Dim<1> >::operator-=(const GridCellIndex<Dim<1> >& rhs) {
  this->mx -= rhs.mx;
  return *this;
}

template<>
inline
GridCellIndex<Dim<2> >&
GridCellIndex<Dim<2> >::operator-=(const GridCellIndex<Dim<2> >& rhs) {
  this->mx -= rhs.mx;
  this->my -= rhs.my;
  return *this;
}

template<>
inline
GridCellIndex<Dim<3> >&
GridCellIndex<Dim<3> >::operator-=(const GridCellIndex<Dim<3> >& rhs) {
  this->mx -= rhs.mx;
  this->my -= rhs.my;
  this->mz -= rhs.mz;
  return *this;
}

//------------------------------------------------------------------------------
// Add an integer to the gridcell indices.
//------------------------------------------------------------------------------
template<>
inline
GridCellIndex<Dim<1> >
GridCellIndex<Dim<1> >::operator+(const int rhs) const {
  return GridCellIndex<Dim<1> >(this->mx + rhs);
}

template<>
inline
GridCellIndex<Dim<2> >
GridCellIndex<Dim<2> >::operator+(const int rhs) const {
  return GridCellIndex<Dim<2> >(this->mx + rhs,
                                this->my + rhs);
}

template<>
inline
GridCellIndex<Dim<3> >
GridCellIndex<Dim<3> >::operator+(const int rhs) const {
  return GridCellIndex<Dim<3> >(this->mx + rhs,
                                this->my + rhs,
                                this->mz + rhs);
}

//------------------------------------------------------------------------------
// Subtract an integer from the gridcell indices.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>
GridCellIndex<Dimension>::operator-(const int /*rhs*/) const {
  VERIFY(false);
  return *this;
}

template<>
inline
GridCellIndex<Dim<1> >
GridCellIndex<Dim<1> >::operator-(const int rhs) const {
  return GridCellIndex<Dim<1> >(this->mx - rhs);
}

template<>
inline
GridCellIndex<Dim<2> >
GridCellIndex<Dim<2> >::operator-(const int rhs) const {
  return GridCellIndex<Dim<2> >(this->mx - rhs,
                                this->my - rhs);
}

template<>
inline
GridCellIndex<Dim<3> >
GridCellIndex<Dim<3> >::operator-(const int rhs) const {
  return GridCellIndex<Dim<3> >(this->mx - rhs,
                                this->my - rhs,
                                this->mz - rhs);
}

//------------------------------------------------------------------------------
// Inline addition with an int.
//------------------------------------------------------------------------------
template<>
inline
GridCellIndex<Dim<1> >&
GridCellIndex<Dim<1> >::operator+=(const int rhs) {
  this->mx += rhs;
  return *this;
}

template<>
inline
GridCellIndex<Dim<2> >&
GridCellIndex<Dim<2> >::operator+=(const int rhs) {
  this->mx += rhs;
  this->my += rhs;
  return *this;
}

template<>
inline
GridCellIndex<Dim<3> >&
GridCellIndex<Dim<3> >::operator+=(const int rhs) {
  this->mx += rhs;
  this->my += rhs;
  this->mz += rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Inline subtraction with an int.
//------------------------------------------------------------------------------
template<>
inline
GridCellIndex<Dim<1> >&
GridCellIndex<Dim<1> >::operator-=(const int rhs) {
  this->mx -= rhs;
  return *this;
}

template<>
inline
GridCellIndex<Dim<2> >&
GridCellIndex<Dim<2> >::operator-=(const int rhs) {
  this->mx -= rhs;
  this->my -= rhs;
  return *this;
}

template<>
inline
GridCellIndex<Dim<3> >&
GridCellIndex<Dim<3> >::operator-=(const int rhs) {
  this->mx -= rhs;
  this->my -= rhs;
  this->mz -= rhs;
  return *this;
}

//------------------------------------------------------------------------------
// Multiply the gridcell indices by an integer.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>
GridCellIndex<Dimension>::operator*(int /*rhs*/) const {
  VERIFY(false);
  return *this;
}

template<>
inline
GridCellIndex<Dim<1> >
GridCellIndex<Dim<1> >::operator*(int rhs) const {
  return GridCellIndex<Dim<1> >(this->mx * rhs);
}

template<>
inline
GridCellIndex<Dim<2> >
GridCellIndex<Dim<2> >::operator*(int rhs) const {
  return GridCellIndex<Dim<2> >(this->mx * rhs,
                                this->my * rhs);
}

template<>
inline
GridCellIndex<Dim<3> >
GridCellIndex<Dim<3> >::operator*(int rhs) const {
  return GridCellIndex<Dim<3> >(this->mx * rhs,
                                this->my * rhs,
                                this->mz * rhs);
}

//------------------------------------------------------------------------------
// Divide the gridcell indices by an integer.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>
GridCellIndex<Dimension>::operator/(int /*rhs*/) const {
  VERIFY(false);
  return *this;
}

template<>
inline
GridCellIndex<Dim<1> >
GridCellIndex<Dim<1> >::operator/(int rhs) const {
  REQUIRE(rhs != 0);
  return GridCellIndex<Dim<1> >(this->mx / rhs - ((this->mx % rhs) < 0 ? 1 : 0));
}

template<>
inline
GridCellIndex<Dim<2> >
GridCellIndex<Dim<2> >::operator/(int rhs) const {
  REQUIRE(rhs != 0);
  return GridCellIndex<Dim<2> >(this->mx / rhs - ((this->mx % rhs) < 0 ? 1 : 0),
                                this->my / rhs - ((this->my % rhs) < 0 ? 1 : 0));
}

template<>
inline
GridCellIndex<Dim<3> >
GridCellIndex<Dim<3> >::operator/(int rhs) const {
  REQUIRE(rhs != 0);
  return GridCellIndex<Dim<3> >(this->mx / rhs - ((this->mx % rhs) < 0 ? 1 : 0),
                                this->my / rhs - ((this->my % rhs) < 0 ? 1 : 0),
                                this->mz / rhs - ((this->mz % rhs) < 0 ? 1 : 0));
}

//------------------------------------------------------------------------------
// Dot product operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::dot(const GridCellIndex<Dimension>& /*rhs*/) const {
  VERIFY(false);
  return 0;
}

template<>
inline
int
GridCellIndex< Dim<1> >::dot(const GridCellIndex< Dim<1> >& rhs) const {
  return this->mx*rhs.mx;
}

template<>
inline
int
GridCellIndex< Dim<2> >::dot(const GridCellIndex< Dim<2> >& rhs) const {
  return this->mx*rhs.mx + this->my*rhs.my;
}

template<>
inline
int
GridCellIndex< Dim<3> >::dot(const GridCellIndex< Dim<3> >& rhs) const {
  return this->mx*rhs.mx + this->my*rhs.my + this->mz*rhs.mz;
}

//------------------------------------------------------------------------------
// Return (-1, 0, 1) if this vector is (less than, equal to, greater than)
// the given vector.
//------------------------------------------------------------------------------
template<>
inline
int
GridCellIndex< Dim<1> >::compare(const GridCellIndex< Dim<1> >& rhs) const {
  return (this->mx < rhs.mx ? -1 :
          this->mx > rhs.mx ?  1 :
          0);
}

template<>
inline
int
GridCellIndex< Dim<2> >::compare(const GridCellIndex< Dim<2> >& rhs) const {
  return (this->my < rhs.my ? -1 :
          this->my > rhs.my ?  1 :
          this->mx < rhs.mx ? -1 :
          this->mx > rhs.mx ?  1 :
          0);
}

template<>
inline
int
GridCellIndex< Dim<3> >::compare(const GridCellIndex< Dim<3> >& rhs) const {
  return (this->mz < rhs.mz ? -1 :
          this->mz > rhs.mz ?  1 :
          this->my < rhs.my ? -1 :
          this->my > rhs.my ?  1 :
          this->mx < rhs.mx ? -1 :
          this->mx > rhs.mx ?  1 :
          0);
}

//------------------------------------------------------------------------------
// Operator ==
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GridCellIndex<Dimension>::operator==(const GridCellIndex<Dimension>& /*rhs*/) const {
  VERIFY(false);
  return false;
}

template<>
inline
bool
GridCellIndex<Dim<1> >::operator==(const GridCellIndex<Dim<1> >& rhs) const {
  return this->mx == rhs.mx;
}

template<>
inline
bool
GridCellIndex<Dim<2> >::operator==(const GridCellIndex<Dim<2> >& rhs) const {
  return ((this->mx == rhs.mx) and
          (this->my == rhs.my));
}

template<>
inline
bool
GridCellIndex<Dim<3> >::operator==(const GridCellIndex<Dim<3> >& rhs) const {
  return ((this->mx == rhs.mx) and
          (this->my == rhs.my) and
          (this->mz == rhs.mz));
}

//------------------------------------------------------------------------------
// Operator !=
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GridCellIndex<Dimension>::operator!=(const GridCellIndex<Dimension>& rhs) const {
  return !(*this == rhs);
}

//------------------------------------------------------------------------------
// Operator <
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GridCellIndex<Dimension>::operator<(const GridCellIndex<Dimension>& rhs) const {
  return this->compare(rhs) == -1;
}

//------------------------------------------------------------------------------
// Operator >
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GridCellIndex<Dimension>::operator>(const GridCellIndex<Dimension>& rhs) const {
  return this->compare(rhs) == 1;
}

//------------------------------------------------------------------------------
// Operator <=
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GridCellIndex<Dimension>::operator<=(const GridCellIndex<Dimension>& rhs) const {
  return this->compare(rhs) <= 0;
}

//------------------------------------------------------------------------------
// Operator >=
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GridCellIndex<Dimension>::operator>=(const GridCellIndex<Dimension>& rhs) const {
  return this->compare(rhs) >= 0;
}

//------------------------------------------------------------------------------
// Operator < GridCellPlane
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GridCellIndex<Dimension>::operator<(const GridCellPlane<Dimension>& rhs) const {
  return rhs > *this;
}

//------------------------------------------------------------------------------
// Operator > GridCellPlane
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GridCellIndex<Dimension>::operator>(const GridCellPlane<Dimension>& rhs) const {
  return rhs < *this;
}

//------------------------------------------------------------------------------
// Operator <= GridCellPlane
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GridCellIndex<Dimension>::operator<=(const GridCellPlane<Dimension>& rhs) const {
  return rhs >= *this;
}

//------------------------------------------------------------------------------
// Operator >= GridCellPlane
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
GridCellIndex<Dimension>::operator>=(const GridCellPlane<Dimension>& rhs) const {
  return rhs <= *this;
}

//------------------------------------------------------------------------------
// Test if the grid cell is contained in the given range.
//------------------------------------------------------------------------------
template<>
inline
bool
GridCellIndex<Dim<1> >::inRange(const GridCellIndex<Dim<1> >& minGridCell,
                                const GridCellIndex<Dim<1> >& maxGridCell) const {
  return this->mx >= minGridCell.mx and this->mx <= maxGridCell.mx;
}

template<>
inline
bool
GridCellIndex<Dim<2> >::inRange(const GridCellIndex<Dim<2> >& minGridCell,
                                const GridCellIndex<Dim<2> >& maxGridCell) const {
  if (this->mx < minGridCell.mx) return false;
  if (this->mx > maxGridCell.mx) return false;

  if (this->my < minGridCell.my) return false;
  if (this->my > maxGridCell.my) return false;

  return true;
}

template<>
inline
bool
GridCellIndex<Dim<3> >::inRange(const GridCellIndex<Dim<3> >& minGridCell,
                                const GridCellIndex<Dim<3> >& maxGridCell) const {
  if (this->mx < minGridCell.mx) return false;
  if (this->mx > maxGridCell.mx) return false;

  if (this->my < minGridCell.my) return false;
  if (this->my > maxGridCell.my) return false;

  if (this->mz < minGridCell.mz) return false;
  if (this->mz > maxGridCell.mz) return false;

  return true;
}

//------------------------------------------------------------------------------
// Return the magnitude as a double.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
GridCellIndex<Dimension>::magnitude() const {
  return std::sqrt((double) this->magnitude2());
}

//------------------------------------------------------------------------------
// Return the square of the magnitude.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::magnitude2() const {
  VERIFY(false);
  return 0;
}

template<>
inline
int
GridCellIndex< Dim<1> >::magnitude2() const {
  return this->mx*this->mx;
}

template<>
inline
int
GridCellIndex< Dim<2> >::magnitude2() const {
  return this->mx*this->mx + this->my*this->my;
}

template<>
inline
int
GridCellIndex< Dim<3> >::magnitude2() const {
  return this->mx*this->mx + this->my*this->my + this->mz*this->mz;
}

//------------------------------------------------------------------------------
// The minimum element.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::minElement() const {
  VERIFY(false);
  return 0;
}

template<>
inline
int
GridCellIndex<Dim<1> >::minElement() const {
  return this->mx;
}

template<>
inline
int
GridCellIndex<Dim<2> >::minElement() const {
  return std::min(this->mx, this->my);
}

template<>
inline
int
GridCellIndex<Dim<3> >::minElement() const {
  return std::min(this->mx, std::min(this->my, this->mz));
}

//------------------------------------------------------------------------------
// The maximum element.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::maxElement() const {
  VERIFY(false);
  return 0;
}

template<>
inline
int
GridCellIndex<Dim<1> >::maxElement() const {
  return this->mx;
}

template<>
inline
int
GridCellIndex<Dim<2> >::maxElement() const {
  return std::max(this->mx, this->my);
}

template<>
inline
int
GridCellIndex<Dim<3> >::maxElement() const {
  return std::max(this->mx, std::max(this->my, this->mz));
}


//------------------------------------------------------------------------------
// The sum of the elements.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::sumElements() const {
  VERIFY(false);
  return 0;
}

template<>
inline
int
GridCellIndex< Dim<1> >::sumElements() const {
  return this->mx;
}

template<>
inline
int
GridCellIndex< Dim<2> >::sumElements() const {
  return this->mx + this->my;
}

template<>
inline
int
GridCellIndex< Dim<3> >::sumElements() const {
  return this->mx + this->my + this->mz;
}

//------------------------------------------------------------------------------
// The product of the elements.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::productElements() const {
  VERIFY(false);
  return 0;
}

template<>
inline
int
GridCellIndex< Dim<1> >::productElements() const {
  return this->mx;
}

template<>
inline
int
GridCellIndex< Dim<2> >::productElements() const {
  return this->mx * this->my;
}

template<>
inline
int
GridCellIndex< Dim<3> >::productElements() const {
  return this->mx * this->my * this->mz;
}

//------------------------------------------------------------------------------
// Return the minimum allowed index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::indexMin() const {
  return mIndexMin;
}

//------------------------------------------------------------------------------
// Return the maximum allowed index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
GridCellIndex<Dimension>::indexMax() const {
  return mIndexMax;
}

//------------------------------------------------------------------------------
// Provide begin/end iterators over the elements of the Vector.
//------------------------------------------------------------------------------
// Non-const versions.
template<typename Dimension>
inline
typename GridCellIndex<Dimension>::iterator
GridCellIndex<Dimension>::begin() {
  return &(this->mx);
}

template<typename Dimension>
inline
typename GridCellIndex<Dimension>::iterator
GridCellIndex<Dimension>::end() {
  return &(this->mx) + Dimension::nDim;
}

// Const versions.
template<typename Dimension>
inline
typename GridCellIndex<Dimension>::const_iterator
GridCellIndex<Dimension>::begin() const {
  return &(this->mx);
}

template<typename Dimension>
inline
typename GridCellIndex<Dimension>::const_iterator
GridCellIndex<Dimension>::end() const {
  return &(this->mx) + Dimension::nDim;
}

//------------------------------------------------------------------------------
// Upack the individual indices into a single int integer value.
//------------------------------------------------------------------------------
// template<>
// inline
// unsigned long long
// GridCellIndex<Dim<1> >::packIndices() const {
//   CHECK(this->mx >= mIndexMin and this->mx < mIndexMax);
//   unsigned long long result = (unsigned long long)(this->mx - mIndexMin);
//   return result;
// }

// template<>
// inline
// unsigned long long
// GridCellIndex<Dim<2> >::packIndices() const {
//   CHECK(this->mx >= mIndexMin and this->mx < mIndexMax);
//   CHECK(this->my >= mIndexMin and this->my < mIndexMax);
//   unsigned long long result = 
//     (unsigned long long)(this->mx - mIndexMin) + 
//     (unsigned long long)(this->my - mIndexMin)*mXMultiplier;
//   return result;
// }

// template<>
// inline
// unsigned long long
// GridCellIndex<Dim<3> >::packIndices() const {
//   CHECK(this->mx >= mIndexMin and this->mx < mIndexMax);
//   CHECK(this->my >= mIndexMin and this->my < mIndexMax);
//   CHECK(this->mz >= mIndexMin and this->mz < mIndexMax);
//   unsigned long long result =
//     (unsigned long long)(this->mx - mIndexMin) + 
//     (unsigned long long)(this->my - mIndexMin)*mXMultiplier +
//     (unsigned long long)(this->mz - mIndexMin)*mXYMultiplier;
//   return result;
// }

//------------------------------------------------------------------------------
// Provide a method to add an integer to the grid cell indices.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>
operator+(int lhs, const GridCellIndex<Dimension>& rhs) {
  GridCellIndex<Dimension> result(rhs);
  return result + lhs;
}

//------------------------------------------------------------------------------
// Provide a method to subtract a grid cell index from an integer.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>
operator-(int lhs, const GridCellIndex<Dimension>& rhs) {
  GridCellIndex<Dimension> result(-rhs);
  return result + lhs;
}

//------------------------------------------------------------------------------
// Provide a method to multiply an integer by a grid cell index.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
GridCellIndex<Dimension>
operator*(int lhs, const GridCellIndex<Dimension>& rhs) {
  GridCellIndex<Dimension> result(rhs);
  return result * lhs;
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
std::ostream&
operator<<(std::ostream& os, const GridCellIndex<Dimension>& gc) {
  os << "GridCell" << Dimension::nDim << "d( ";
  for (int i = 0; i != Dimension::nDim; ++i) os << gc(i) << " ";
  os << ")";
  return os;
}

}
