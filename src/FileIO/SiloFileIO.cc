//---------------------------------Spheral++----------------------------------//
// SiloFileIO -- Provide the interface to silo file objects.
//
// Created by JMO, Sat Feb  7 23:06:03 PST 2015
//----------------------------------------------------------------------------//
#include "SiloFileIO.hh"
#include "Field/Field.hh"

#include "Utilities/DBC.hh"
#include "boost/algorithm/string.hpp"
#include "boost/algorithm/string/replace.hpp"

#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

namespace {
  
//------------------------------------------------------------------------------
// Mangle a path name for silo consumption.
//------------------------------------------------------------------------------
string SILO_mangle(const string& x) {
  string result = x;
  boost::replace_all(result, " ", "_");
  boost::replace_all(result, ".", "_p_");
  boost::replace_all(result, ",", "_c_");
  boost::replace_all(result, "|", "_P_");
  boost::replace_all(result, "-", "_d_");
  return result;
}

//------------------------------------------------------------------------------
// Set the directory in the current silo file to the path portion of a path,
// and return the variable name section.
//------------------------------------------------------------------------------
string setdir(DBfile* filePtr, const string& ipath) {

  // Mangle the path into something acceptable to silo.
  string path = SILO_mangle(ipath);

  // We always start from the top.
  VERIFY2(DBSetDir(filePtr, "/") == 0,
          "SiloFileIO ERROR: unable to change to path /");

  // Split up the path into directories and the var name.
  vector<string> components;
  boost::split(components, path, boost::is_any_of("/"));
  CHECK(components.size() > 0);

  // Set the path and return just the variable name.
  for (unsigned i = 0; i != components.size() - 1; ++i) {
    const string dirName = components[i];
    if (dirName.size() > 0) {
      // Check if this directory already exists or not.
      DBtoc* toc = DBGetToc(filePtr);
      bool exists = false;
      int j = 0;
      while (not exists and j != toc->ndir) {
        exists = (dirName == string(toc->dir_names[j]));
        ++j;
      }
      if (not exists) {
        // std::cerr << " --> " << path << std::endl
        //           << " --> " << dirName.size() << std::endl
        //           << " --> " << dirName << std::endl;
        VERIFY2(DBMkDir(filePtr, dirName.c_str()) == 0,
                "SiloFileIO ERROR: unable to create path " << dirName << " of " << path);
      }
      VERIFY2(DBSetDir(filePtr, dirName.c_str()) == 0,
              "SiloFileIO ERROR: unable to change to path " << dirName << " of " << path);
    }
  }

  return components.back();
}


//------------------------------------------------------------------------------
// Common methods for reading/writing a lowly int to a silo file.
//------------------------------------------------------------------------------
void writeInt(DBfile* filePtr, const int& value, const string path) {
  const string varname = setdir(filePtr, path);
  int dims[1] = {1};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(filePtr, varname.c_str(), &value, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write int variable " << path);
}

void readInt(DBfile* filePtr, int& value, const string path) {
  // const string varname = setdir(mFilePtr, path);
  // CHECK2(DBReadVar(filePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << path);
  const string varname = setdir(filePtr, path);
  value = *static_cast<int*>(DBGetVar(filePtr, varname.c_str()));
}

//------------------------------------------------------------------------------
// Common methods for reading/writing a std::string to a silo file.
//------------------------------------------------------------------------------
void writeString(DBfile* filePtr, const string& value, const string path) {
  const int size = value.size();
  int dims[1] = {size};
  writeInt(filePtr, dims[0], path + "/size");
  if (dims[0] > 0) {
    const string varname = setdir(filePtr, path + "/value");
    VERIFY2(DBWrite(filePtr, varname.c_str(), value.c_str(), dims, 1, DB_CHAR) == 0,
            "SiloFileIO ERROR: unable to write string variable " << path);
  }
}

void readString(DBfile* filePtr, string& value, const string path) {
  int valsize;
  readInt(filePtr, valsize, path + "/size");
  if (valsize == 0) {
    value = "";
  } else {
    std::vector<char> cvalue(valsize + 1);  // Do we need to allow space for trailing null?
    const string varname = setdir(filePtr, path + "/value");
    VERIFY2(DBReadVar(filePtr, varname.c_str(), &cvalue[0]) == 0,
            "SiloFileIO ERROR: failed to read string variable " << path);
    value = string(&cvalue[0], &cvalue[valsize]);
  }
}

}   // anonymous namespace

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
SiloFileIO::SiloFileIO():
  FileIO(),
  mFilePtr(0) {
}

//------------------------------------------------------------------------------
// Construct and open the given file.
//------------------------------------------------------------------------------
SiloFileIO::
SiloFileIO(const string fileName, AccessType access):
  FileIO(fileName, access),
  mFilePtr(0) {
  open(fileName, access);
  ENSURE(mFileOpen && mFilePtr != 0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
SiloFileIO::~SiloFileIO() {
  close();
}

//------------------------------------------------------------------------------
// Open a SiloFile file with the specified access.
//------------------------------------------------------------------------------
void
SiloFileIO::open(const string fileName, AccessType access) {
  VERIFY2(mFilePtr == 0 and mFileOpen == false,
          "ERROR: attempt to reopen SiloFileIO object.");

  string fullFileName = fileName;
  if (fullFileName.find(".silo") == string::npos) {
    fullFileName += ".silo";
  }

  if (access == AccessType::Read) {
    mFilePtr = DBOpen(fullFileName.c_str(), DB_HDF5, DB_READ);
  } else {
    mFilePtr = DBCreate(fullFileName.c_str(), DB_CLOBBER, DB_LOCAL, "Spheral++ restart file.", DB_HDF5);
  }
  VERIFY2(mFilePtr != 0, "SiloFileIO ERROR: unable to open " << fullFileName);
  mFileOpen = true;
}

//------------------------------------------------------------------------------
// Close the current file.
//------------------------------------------------------------------------------
void
SiloFileIO::close() {
  if (mFilePtr != 0) {
    VERIFY2(DBClose(mFilePtr) == 0,
            "SiloFileIO ERROR: unable to close file.");
    mFilePtr = 0;
  }
  mFileOpen = false;
}

//------------------------------------------------------------------------------
// Check if the specified path is in the file.
//------------------------------------------------------------------------------
bool
SiloFileIO::pathExists(const std::string path) const {
  const auto varname = setdir(mFilePtr, path);
  return DBInqVarExists(mFilePtr, varname.c_str()) == 1;
}

//------------------------------------------------------------------------------
// Write an unsigned to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const unsigned& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {1};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a size_t to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const size_t& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {1};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const int& value, const string path) {
  writeInt(mFilePtr, value, path);
}

//------------------------------------------------------------------------------
// Write a bool to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const bool& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {1};
  int ivalue = value ? 1 : 0;
  CONTRACT_VAR(dims);
  CONTRACT_VAR(ivalue);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &ivalue, dims, 1, DB_INT) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a double to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const double& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {1};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &value, dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a string to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const string& value, const string path) {
  writeString(mFilePtr, value, path);
}

//------------------------------------------------------------------------------
// Write a vector<int> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const std::vector<int>& value, const string path) {
  const int size = value.size();
  writeInt(mFilePtr, size, path + "/size");
  if (size > 0) {
    const string varname = setdir(mFilePtr, path + "/value");
    int dims[1] = {size};
    CONTRACT_VAR(dims);
    VERIFY2(DBWrite(mFilePtr, varname.c_str(), static_cast<void*>(const_cast<int*>(&value[0])), dims, 1, DB_INT) == 0,
            "SiloFileIO ERROR: unable to write std::vector " << path);
  }
}

//------------------------------------------------------------------------------
// Write a vector<double> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const std::vector<double>& value, const string path) {
  const int size = value.size();
  writeInt(mFilePtr, size, path + "/size");
  if (size > 0) {
    const string varname = setdir(mFilePtr, path + "/value");
    int dims[1] = {size};
    CONTRACT_VAR(dims);
    VERIFY2(DBWrite(mFilePtr, varname.c_str(), static_cast<void*>(const_cast<double*>(&value[0])), dims, 1, DB_DOUBLE) == 0,
            "SiloFileIO ERROR: unable to write std::vector " << path);
  }
}

//------------------------------------------------------------------------------
// Write a vector<string> to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const std::vector<string>& value, const string path) {
  const unsigned n = value.size();
  vector<int> dim_stuff(n);
  string stuff;
  for (unsigned i = 0; i != n; ++i) {
    dim_stuff[i] = value[i].size();
    stuff += value[i];
  }
  this->write(dim_stuff, path + "/dim_stuff");
  this->write(stuff, path + "/stuff");
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::Vector& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<1>::Vector::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::Tensor& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<1>::Tensor::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::SymTensor& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<1>::SymTensor::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<1>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<1>::ThirdRankTensor& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<1>::ThirdRankTensor::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::Vector& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<2>::Vector::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::Tensor& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<2>::Tensor::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::SymTensor& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<2>::SymTensor::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<2>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<2>::ThirdRankTensor& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<2>::ThirdRankTensor::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Vector to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::Vector& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<3>::Vector::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::Tensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::Tensor& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<3>::Tensor::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::SymTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::SymTensor& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<3>::SymTensor::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Write a Dim<3>::ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::write(const Dim<3>::ThirdRankTensor& value, const string path) {
  const string varname = setdir(mFilePtr, path);
  int dims[1] = {int(Dim<3>::ThirdRankTensor::numElements())};
  CONTRACT_VAR(dims);
  VERIFY2(DBWrite(mFilePtr, varname.c_str(), &(*value.begin()), dims, 1, DB_DOUBLE) == 0,
          "SiloFileIO ERROR: unable to write variable " << path);
}

//------------------------------------------------------------------------------
// Read an unsigned from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(unsigned& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << path);
  value = *static_cast<unsigned*>(DBGetVar(mFilePtr, varname.c_str()));
}

//------------------------------------------------------------------------------
// Read a size_t from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(size_t& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << path);
  value = *static_cast<size_t*>(DBGetVar(mFilePtr, varname.c_str()));
}

//------------------------------------------------------------------------------
// Read an int to the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(int& value, const string path) const {
  readInt(mFilePtr, value, path);
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(bool& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  const int ivalue = *static_cast<int*>(DBGetVar(mFilePtr, varname.c_str()));
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &ivalue) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << path);
  value = (ivalue == 1 ? true : false);
}

//------------------------------------------------------------------------------
// Read a double from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(double& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  value = *static_cast<double*>(DBGetVar(mFilePtr, varname.c_str()));
  // VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &value) == 0,
  //         "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a string from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(string& value, const string path) const {
  readString(mFilePtr, value, path);
}

//------------------------------------------------------------------------------
// Read a vector<int> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(std::vector<int>& value, const string path) const {
  int size;
  readInt(mFilePtr, size, path + "/size");
  value.resize(size);
  if (size > 0) {
    const string varname = setdir(mFilePtr, path + "/value");
    VERIFY2(DBReadVar(mFilePtr, varname.c_str(), static_cast<void*>(&value[0])) == 0,
            "SiloFileIO ERROR: unable to read std::vector " << path);
  }
}

//------------------------------------------------------------------------------
// Read a vector<double> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(std::vector<double>& value, const string path) const {
  int size;
  readInt(mFilePtr, size, path + "/size");
  value.resize(size);
  if (size > 0) {
    const string varname = setdir(mFilePtr, path + "/value");
    VERIFY2(DBReadVar(mFilePtr, varname.c_str(), static_cast<void*>(&value[0])) == 0,
            "SiloFileIO ERROR: unable to read std::vector " << path);
  }
}

//------------------------------------------------------------------------------
// Read a vector<string> from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(vector<string>& value, const string path) const {
  vector<int> dim_stuff;
  string stuff;
  this->read(dim_stuff, path + "/dim_stuff");
  this->read(stuff, path + "/stuff");
  const unsigned n = dim_stuff.size();
  value = vector<string>(n);
  unsigned i = 0;
  for (unsigned k = 0; k != n; ++k) {
    value[k] = stuff.substr(i, dim_stuff[k]);
    i += dim_stuff[k];
  }
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::Vector& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::Tensor& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::SymTensor& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<1>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<1>::ThirdRankTensor& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::Vector& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::Tensor& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::SymTensor& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<2>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<2>::ThirdRankTensor& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Vector from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::Vector& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::Tensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::Tensor& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::SymTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::SymTensor& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

//------------------------------------------------------------------------------
// Read a Dim<3>::ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
SiloFileIO::read(Dim<3>::ThirdRankTensor& value, const string path) const {
  const string varname = setdir(mFilePtr, path);
  VERIFY2(DBReadVar(mFilePtr, varname.c_str(), &(*value.begin())) == 0,
          "SiloFileIO ERROR: unable to read variable " << path);
}

}
