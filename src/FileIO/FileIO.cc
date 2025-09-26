//---------------------------------Spheral++----------------------------------//
// FileIO -- Provide the generic interface to file objects.
//
// Created by J. Michael Owen, Wed Feb  7 22:59:14 PST 2001
//----------------------------------------------------------------------------//
#include "FileIO.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "config.hh"

#include "boost/algorithm/string.hpp"

#include <cstring>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

#ifdef SPHERAL_ENABLE_PYTHON
namespace py = pybind11;
#endif

namespace Spheral {

//------------------------------------------------------------------------------
// Empty Constructor.
//------------------------------------------------------------------------------
FileIO::FileIO():
  mFileName(""),
  mAccess(AccessType::Undefined),
  mFileOpen(false) {
}

//------------------------------------------------------------------------------
// Construct with the given file name and access type.
//------------------------------------------------------------------------------
FileIO::FileIO(const string filename, AccessType access):
  mFileName(filename),
  mAccess(access),
  mFileOpen(false) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
FileIO::~FileIO() {
}

//------------------------------------------------------------------------------
// Split the given path into its component strings.
//------------------------------------------------------------------------------
vector<string>
FileIO::splitPathComponents(const string path) const {
  vector<string> components;
  boost::split(components, path, boost::is_any_of("/"));
  return components;
}

// Reverse operation
std::string
FileIO::joinPathComponents(const std::vector<std::string>& components) const {
  string result = "";
  for (const auto& s: components) result += "/" + s;
  return result;
}

//------------------------------------------------------------------------------
// Write GeomPlanes.
//------------------------------------------------------------------------------
void
FileIO::write(const GeomPlane<Dim<1> >& value, const string path) {
  write(value.point(), path + "/point");
  write(value.normal(), path + "/normal");
}

void
FileIO::write(const GeomPlane<Dim<2> >& value, const string path) {
  write(value.point(), path + "/point");
  write(value.normal(), path + "/normal");
}

void
FileIO::write(const GeomPlane<Dim<3> >& value, const string path) {
  write(value.point(), path + "/point");
  write(value.normal(), path + "/normal");
}

//------------------------------------------------------------------------------
// Read GeomPlanes.
//------------------------------------------------------------------------------
void
FileIO::read(GeomPlane<Dim<1> >& value, const string path) const {
  Dim<1>::Vector point, normal;
  read(point, path + "/point");
  read(normal, path + "/normal");
  value.point(point);
  value.normal(normal);
}

void
FileIO::read(GeomPlane<Dim<2> >& value, const string path) const {
  Dim<2>::Vector point, normal;
  read(point, path + "/point");
  read(normal, path + "/normal");
  value.point(point);
  value.normal(normal);
}

void
FileIO::read(GeomPlane<Dim<3> >& value, const string path) const {
  Dim<3>::Vector point, normal;
  read(point, path + "/point");
  read(normal, path + "/normal");
  value.point(point);
  value.normal(normal);
}

//------------------------------------------------------------------------------
// Write polytopes
//------------------------------------------------------------------------------
void
FileIO::write(const Dim<1>::FacetedVolume& value, const string path) {
  std::vector<char> buf;
  packElement(value, buf);
  this->write_vector_char(buf, path);
}

void
FileIO::write(const Dim<2>::FacetedVolume& value, const string path) {
  std::vector<char> buf;
  packElement(value, buf);
  this->write_vector_char(buf, path);
}

void
FileIO::write(const Dim<3>::FacetedVolume& value, const string path) {
  std::vector<char> buf;
  packElement(value, buf);
  this->write_vector_char(buf, path);
}

//------------------------------------------------------------------------------
// Read polytopes
//------------------------------------------------------------------------------
void
FileIO::read(Dim<1>::FacetedVolume& value, const string path) const {
  const auto buf = this->read_vector_char(path);
  auto itr = buf.begin();
  unpackElement(value, itr, buf.end());
  ENSURE(itr == buf.end());
}

void
FileIO::read(Dim<2>::FacetedVolume& value, const string path) const {
  const auto buf = this->read_vector_char(path);
  auto itr = buf.begin();
  unpackElement(value, itr, buf.end());
  ENSURE(itr == buf.end());
}

void
FileIO::read(Dim<3>::FacetedVolume& value, const string path) const {
  const auto buf = this->read_vector_char(path);
  auto itr = buf.begin();
  unpackElement(value, itr, buf.end());
  ENSURE(itr == buf.end());
}

//------------------------------------------------------------------------------
// Write uniform_random
//------------------------------------------------------------------------------
void
FileIO::write(const uniform_random& value, const string path) {
  write((unsigned) value.seed(), path + "/seed");
  write((unsigned) value.numCalls(), path + "/numCalls");
  write(value.min(), path + "/min");
  write(value.max(), path + "/max");
}

//------------------------------------------------------------------------------
// Read uniform_random
//------------------------------------------------------------------------------
void
FileIO::read(uniform_random& value, const string path) const {
  unsigned seed, numCalls;
  double a, b;
  read(seed, path + "/seed");
  read(numCalls, path + "/numCalls");
  read(a, path + "/min");
  read(b, path + "/max");
  value.seed(seed);
  value.range(a, b);
  value.advance(numCalls);
}

//------------------------------------------------------------------------------
// Provide access to the string write method with char*
//------------------------------------------------------------------------------
void
FileIO::write(const char* value, const string path) {
  write(string(value, strlen(value)), path);
}

//------------------------------------------------------------------------------
// Provide access to the string read method with char*
//------------------------------------------------------------------------------
void
FileIO::read(char* value, const string path) const {
  string strValue;
  read(strValue, path);
  strcpy(value, strValue.c_str());
}

//------------------------------------------------------------------------------
// Return the group (directory) component of the given path.
//------------------------------------------------------------------------------
string
FileIO::groupName(const string path) const {

  const vector<string> components = splitPathComponents(path);
  CHECK(components.size() > 0);

  string groupName = "";
  for (vector<string>::const_iterator nameItr = components.begin();
       nameItr < components.end() - 1;
       ++nameItr) groupName += *nameItr;

  return groupName;
}

//------------------------------------------------------------------------------
// Return the variable name component of the given path.
//------------------------------------------------------------------------------
string
FileIO::variableName(const string path) const {
  const vector<string> components = splitPathComponents(path);
  CHECK(components.size() > 0);
  return components[components.size() - 1];
}

//------------------------------------------------------------------------------
// Return the current file name.
//------------------------------------------------------------------------------
const string&
FileIO::fileName() const {
  return mFileName;
}

//------------------------------------------------------------------------------
// Return the access being specified.
//------------------------------------------------------------------------------
AccessType
FileIO::access() const {
  return mAccess;
}

//------------------------------------------------------------------------------
// Return the flag indicating whether the file is open or not.
//------------------------------------------------------------------------------
bool
FileIO::fileOpen() const {
  return mFileOpen;
}

#ifdef SPHERAL_ENABLE_PYTHON
//------------------------------------------------------------------------------
// Python objects (handle with pickle)
//------------------------------------------------------------------------------
void
FileIO::writeObject(py::object thing, const std::string path) {
  auto pickle = py::module_::import("pickle");
  auto dumps = pickle.attr("dumps");
  auto stuff = dumps(thing);
  this->writeBytes(stuff, path);
}

py::object
FileIO::readObject(const std::string path) const {
  auto pickle = py::module_::import("pickle");
  auto loads = pickle.attr("loads");
  auto stuff = this->readBytes(path);
  return loads(stuff);
}

//------------------------------------------------------------------------------
// Python bytes
//------------------------------------------------------------------------------
void
FileIO::writeBytes(py::bytes stuff, const std::string path) {
  std::string stuff_str = std::string(stuff);
  std::vector<char> buf(stuff_str.begin(), stuff_str.end());
  ENSURE(buf.size() == py::len(stuff));
  this->write_vector_char(buf, path);
}

py::bytes
FileIO::readBytes(const std::string path) const {
  auto buf = this->read_vector_char(path);
  return py::bytes(std::string(buf.begin(), buf.end()));
}
#endif

}
