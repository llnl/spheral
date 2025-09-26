//---------------------------------Spheral++----------------------------------//
// FileIO -- Provide the generic interface to file objects.
//
// Created by JMO, Tue Jul 11 22:19:37 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral__FileIO_hh__
#define __Spheral__FileIO_hh__

#include "Geometry/Dimension.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/uniform_random.hh"
#include "Utilities/packElement.hh"

#include <vector>
#include <string>
#include <sstream>
#include "config.hh"

#ifdef SPHERAL_ENABLE_PYTHON
#include "pybind11/pybind11.h"
#include "Utilities/SPHERAL_DLL_EXPORT.hh"
#endif

namespace Spheral {

template<typename Dimension> class GeomPlane;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> struct RKCoefficients;

// Define the standard file access types.
enum class AccessType {
  Undefined = -1,
  Create = 0,
  Read = 1,
  Write = 2,
  ReadWrite = 3
};

class FileIO {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  FileIO();
  FileIO(const std::string filename, AccessType access);

  // Destructor.
  virtual ~FileIO();

  // All File objects must provide methods to open and close the files.
  virtual void open(const std::string fileName, AccessType access) = 0;
  virtual void close() = 0;

  //******************************************************************************
  // Methods all FileIO descendent classes must provide.
  //******************************************************************************
  // Check if the specified path is in the file.
  virtual bool pathExists(const std::string path) const = 0;

  // All FileIO objects had better be able to read and write the primitive 
  // DataTypes.
  virtual void write(const unsigned& value, const std::string path) = 0;
  virtual void write(const size_t& value, const std::string path) = 0;
  virtual void write(const int& value, const std::string path) = 0;
  virtual void write(const bool& value, const std::string path) = 0;
  virtual void write(const double& value, const std::string path) = 0;
  virtual void write(const std::string& value, const std::string path) = 0;
  virtual void write(const std::vector<int>& value, const std::string path) = 0;
  virtual void write(const std::vector<double>& value, const std::string path) = 0;
  virtual void write(const std::vector<std::string>& value, const std::string path) = 0;

  virtual void write(const Dim<1>::Vector& value, const std::string path) = 0;
  virtual void write(const Dim<1>::Tensor& value, const std::string path) = 0;
  virtual void write(const Dim<1>::SymTensor& value, const std::string path) = 0;
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string path) = 0;

  virtual void write(const Dim<2>::Vector& value, const std::string path) = 0;
  virtual void write(const Dim<2>::Tensor& value, const std::string path) = 0;
  virtual void write(const Dim<2>::SymTensor& value, const std::string path) = 0;
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string path) = 0;

  virtual void write(const Dim<3>::Vector& value, const std::string path) = 0;
  virtual void write(const Dim<3>::Tensor& value, const std::string path) = 0;
  virtual void write(const Dim<3>::SymTensor& value, const std::string path) = 0;
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string path) = 0;

  virtual void read(unsigned& value, const std::string path) const = 0;
  virtual void read(size_t& value, const std::string path) const = 0;
  virtual void read(int& value, const std::string path) const = 0;
  virtual void read(bool& value, const std::string path) const = 0;
  virtual void read(double& value, const std::string path) const = 0;
  virtual void read(std::string& value, const std::string path) const = 0;
  virtual void read(std::vector<int>& value, const std::string path) const = 0;
  virtual void read(std::vector<double>& value, const std::string path) const = 0;
  virtual void read(std::vector<std::string>& value, const std::string path) const = 0;

  virtual void read(Dim<1>::Vector& value, const std::string path) const = 0;
  virtual void read(Dim<1>::Tensor& value, const std::string path) const = 0;
  virtual void read(Dim<1>::SymTensor& value, const std::string path) const = 0;
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string path) const = 0;

  virtual void read(Dim<2>::Vector& value, const std::string path) const = 0;
  virtual void read(Dim<2>::Tensor& value, const std::string path) const = 0;
  virtual void read(Dim<2>::SymTensor& value, const std::string path) const = 0;
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string path) const = 0;

  virtual void read(Dim<3>::Vector& value, const std::string path) const = 0;
  virtual void read(Dim<3>::Tensor& value, const std::string path) const = 0;
  virtual void read(Dim<3>::SymTensor& value, const std::string path) const = 0;
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string path) const = 0;

  //******************************************************************************
  // These methods are useful for the primitive types that are problematic
  // to return by reference from python.
  virtual void write_unsigned_int(const unsigned value, const std::string path)                                    { this->write(value, path); }
  virtual void write_size_t(const size_t value, const std::string path)                                            { this->write(value, path); }
  virtual void write_int(const int value, const std::string path)                                                  { this->write(value, path); }
  virtual void write_bool(const bool value, const std::string path)                                                { this->write(value, path); }
  virtual void write_double(const double value, const std::string path)                                            { this->write(value, path); }
  virtual void write_string(const std::string value, const std::string path)                                       { this->write(value, path); }
  virtual void write_vector_char(const std::vector<char>& value,
                                 const std::string path)                                                           { this->write(std::string(value.begin(), value.end()), path); }

  virtual unsigned read_unsigned_int(const std::string path) const                                                 { unsigned result;    this->read(result, path); return result; }
  virtual size_t read_size_t(const std::string path) const                                                         { size_t result;      this->read(result, path); return result; }
  virtual int read_int(const std::string path) const                                                               { int result;         this->read(result, path); return result; }
  virtual bool read_bool(const std::string path) const                                                             { bool result;        this->read(result, path); return result; }
  virtual double read_double(const std::string path) const                                                         { double result;      this->read(result, path); return result; }
  virtual std::string read_string(const std::string path) const                                                    { std::string result; this->read(result, path); return result; }
  virtual std::vector<char> read_vector_char(const std::string path) const                                         { auto result = this->read_string(path); return std::vector<char>(result.begin(), result.end()); }

  // Fields
  template<typename Dimension, typename DataType> void write(const Field<Dimension, DataType>& value, const std::string path);
  template<typename Dimension, typename DataType> void read(Field<Dimension, DataType>& value, const std::string path) const;

  // FieldLists
  template<typename Dimension, typename DataType> void write(const FieldList<Dimension, DataType>& value, const std::string path);
  template<typename Dimension, typename DataType> void read(FieldList<Dimension, DataType>& value, const std::string path) const;

  // Read/write planes using the primitive methods.
  void write(const GeomPlane<Dim<1> >& value, const std::string path);
  void write(const GeomPlane<Dim<2> >& value, const std::string path);
  void write(const GeomPlane<Dim<3> >& value, const std::string path);

  void read(GeomPlane<Dim<1> >& value, const std::string path) const;
  void read(GeomPlane<Dim<2> >& value, const std::string path) const;
  void read(GeomPlane<Dim<3> >& value, const std::string path) const;

  // Read/write polytopes
  void write(const Dim<1>::FacetedVolume& value, const std::string path);
  void write(const Dim<2>::FacetedVolume& value, const std::string path);
  void write(const Dim<3>::FacetedVolume& value, const std::string path);

  void read(Dim<1>::FacetedVolume& value, const std::string path) const;
  void read(Dim<2>::FacetedVolume& value, const std::string path) const;
  void read(Dim<3>::FacetedVolume& value, const std::string path) const;

  // Provide char* read/write methods, which will simply call the string
  // methods provided by descendents.
  void write(const char* value, const std::string path);
  void read(char* value, const std::string path) const;

  // Read/write uniform_random
  void write(const uniform_random& value, const std::string path);
  void read(uniform_random& value, const std::string path) const;

  // Write/read a vector<DataType> if DataType is a primitive we already know about.
  template<typename DataType> void write(const std::vector<DataType>& x, const std::string path);
  template<typename DataType> void read(std::vector<DataType>& x, const std::string path) const;

  // Helper functions to split/join a string up into substrings delimited by '/'.
  std::vector<std::string> splitPathComponents(const std::string path) const;
  std::string joinPathComponents(const std::vector<std::string>& components) const;

  // Return the group (directory) component of a path.
  std::string groupName(const std::string path) const;

  // Return the variable component of a path.
  std::string variableName(const std::string path) const;

  // Allow const access to some of the member data.
  const std::string& fileName() const;
  AccessType access() const;
  bool fileOpen() const;

  // Safe method to try and read from a path if it exists.
  // Returns: 0 => successful
  //          1 => path does not exist
  //          2 => unable to read value
  template<typename T> int readIfAvailable(T& value, const std::string path) const;

#ifdef SPHERAL_ENABLE_PYTHON
  // PyObjects for Python
  SPHERAL_DLL_PUBLIC virtual void writeObject(pybind11::object thing, const std::string path);
  SPHERAL_DLL_PUBLIC virtual pybind11::object readObject(const std::string path) const;

  // pybind11::bytes
  SPHERAL_DLL_PUBLIC virtual void writeBytes(pybind11::bytes thing, const std::string path);
  SPHERAL_DLL_PUBLIC virtual pybind11::bytes readBytes(const std::string path) const;

  // When python tries to read/write bytes objects, use the correct method to avoid casting to a string
  SPHERAL_DLL_PUBLIC virtual void write(pybind11::bytes& thing, const std::string path)       { this->writeBytes(thing, path); }
  SPHERAL_DLL_PUBLIC virtual void read(pybind11::bytes& thing, const std::string path)  const { thing = this->readBytes(path); }
#endif

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Descendent class are allowed to directly diddle this common data.
  std::string mFileName;
  AccessType mAccess;
  bool mFileOpen;

private:
  //--------------------------- Private Interface ---------------------------//
  // Don't allow assignment.
  FileIO& operator=(const FileIO& rhs);

  // Private methods to help with std::vector specializations
  template<typename Value> void writeVector(const std::vector<Value>& x, const std::string path);
  template<typename Value> void readVector(std::vector<Value>& x, const std::string path) const;
};

}

#include "FileIOInline.hh"

#endif
