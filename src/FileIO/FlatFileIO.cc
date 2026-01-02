//---------------------------------Spheral++----------------------------------//
// FlatFileIO -- Provide the interface to FlatFile file objects.
//
// Created by JMO, Fri Apr 13 01:19:02 PDT 2001
//----------------------------------------------------------------------------//
#include "FlatFileIO.hh"
#include "Field/Field.hh"

using std::vector;
using std::string;
using std::pair;
using std::make_pair;
using std::ios;

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
FlatFileIO::FlatFileIO():
  FileIO(),
  mPrecision(20),
  mFilePtr(0) {
}

//------------------------------------------------------------------------------
// Construct and open the given file.
//------------------------------------------------------------------------------
FlatFileIO::
FlatFileIO(const string fileName, AccessType access, FlatFileFormat format):
  FileIO(fileName, access),
  mPrecision(20),
  mFilePtr(0),
  mFileFormat(format) {
  open(fileName, access);
  ENSURE(mFileOpen && mFilePtr != 0);
}

// Do the same thing with an int in the place of the enum, to help Pyffle.
FlatFileIO::
FlatFileIO(const string fileName, int access, int format):
  FileIO(fileName, AccessType(access)),
  mPrecision(20),
  mFilePtr(0),
  mFileFormat(FlatFileFormat(format)) {
  open(fileName, AccessType(access));
  ENSURE(mFileOpen && mFilePtr != 0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
FlatFileIO::~FlatFileIO() {
  close();
}

//------------------------------------------------------------------------------
// Open a FlatFile file with the specified access.
//------------------------------------------------------------------------------
void
FlatFileIO::open(const string fileName, AccessType access) {

  // If we currently have a file open and attached to this object, close it!
  close();
  CHECK(mFilePtr == 0 && mFileOpen == false);

  // Build the file opening mode.
  ios::openmode mode;
  switch(access) {
  case AccessType::Create:
    mode = ios::out; //ios::trunc;
    break;
  case AccessType::Read:
    mode = ios::in;
    break;
  case AccessType::Write:
    mode = ios::out;
    break;
  case AccessType::ReadWrite:
    mode = ios::in | ios::out;
    break;
  default:
    VERIFY2(false, "Unhandled case in switch!");
  }

  if (mFileFormat == FlatFileFormat::binary) mode = mode | ios::binary;

  // Open a file stream and attach it to this objects pointer.
  mFilePtr = new std::fstream(fileName.c_str(), mode);
  mFileOpen = mFilePtr->is_open();
  ENSURE(mFilePtr != 0 && mFileOpen == true);

  // Set the precision we will be using to write to the file.
  (*mFilePtr).precision(mPrecision);
}

//------------------------------------------------------------------------------
// Close the current file.
//------------------------------------------------------------------------------
void
FlatFileIO::close() {

  if (mFilePtr != 0) {
    mFilePtr->flush();
    mFilePtr->close();
    delete mFilePtr;
    mFilePtr = 0;
  }
  mFileOpen = false;
}

//------------------------------------------------------------------------------
// Check if the specified path is in the file.
//------------------------------------------------------------------------------
bool
FlatFileIO::pathExists(const std::string path) const {
  try {
    findPathName(path);
    return true;
  } catch(...) {
    return false;
  }
}

//------------------------------------------------------------------------------
// Write an unsigned to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::write(const unsigned& value, const string path) {
  writeGenericType(value, path);
}

//------------------------------------------------------------------------------
// Write a size_t to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::write(const size_t& value, const string path) {
  writeGenericType(value, path);
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::write(const int& value, const string path) {
  writeGenericType(value, path);
}

//------------------------------------------------------------------------------
// Write a bool to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::write(const bool& value, const string path) {
  writeGenericType(value, path);
}

//------------------------------------------------------------------------------
// Write a Scalar to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const double& value, const string path) {
  writeGenericType(value, path);
}

//------------------------------------------------------------------------------
// Write a Vector to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const Dim<1>::Vector& value, const string path) {
  writeGenericType(value, path);
}

void
FlatFileIO::
write(const Dim<2>::Vector& value, const string path) {
  writeGenericType(value, path);
}

void
FlatFileIO::
write(const Dim<3>::Vector& value, const string path) {
  writeGenericType(value, path);
}

//------------------------------------------------------------------------------
// Write a Tensor to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const Dim<1>::Tensor& value, const string path) {
  writeGenericType(value, path);
}

void
FlatFileIO::
write(const Dim<2>::Tensor& value, const string path) {
  writeGenericType(value, path);
}

void
FlatFileIO::
write(const Dim<3>::Tensor& value, const string path) {
  writeGenericType(value, path);
}

//------------------------------------------------------------------------------
// Write a SymTensor to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const Dim<1>::SymTensor& value, const string path) {
  writeGenericType(value, path);
}

void
FlatFileIO::
write(const Dim<2>::SymTensor& value, const string path) {
  writeGenericType(value, path);
}

void
FlatFileIO::
write(const Dim<3>::SymTensor& value, const string path) {
  writeGenericType(value, path);
}

//------------------------------------------------------------------------------
// Write a ThirdRankTensor to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const Dim<1>::ThirdRankTensor& value, const string path) {
  writeGenericType(value, path);
}

void
FlatFileIO::
write(const Dim<2>::ThirdRankTensor& value, const string path) {
  writeGenericType(value, path);
}

void
FlatFileIO::
write(const Dim<3>::ThirdRankTensor& value, const string path) {
  writeGenericType(value, path);
}

//------------------------------------------------------------------------------
// Write a string to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const string& value, const string path) {
  writeGenericType(value, path);
}

//------------------------------------------------------------------------------
// Write a std::vector<int> to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const std::vector<int>& value, const string path) {
  std::vector<char> buf;
  packElement(value, buf);
  std::string strbuf(buf.begin(), buf.end());
  writeGenericType(strbuf, path);
}

//------------------------------------------------------------------------------
// Write a std::vector<double> to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const std::vector<double>& value, const string path) {
  std::vector<char> buf;
  packElement(value, buf);
  std::string strbuf(buf.begin(), buf.end());
  writeGenericType(strbuf, path);
}

//------------------------------------------------------------------------------
// Write a std::vector<string> to the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
write(const std::vector<string>& value, const string path) {
  std::vector<char> buf;
  packElement(value, buf);
  std::string strbuf(buf.begin(), buf.end());
  writeGenericType(strbuf, path);
}

//------------------------------------------------------------------------------
// Read an unsigned from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(unsigned& value, const string path) const {
  readGenericType(value, path);
}

//------------------------------------------------------------------------------
// Read a size_t from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(size_t& value, const string path) const {
  readGenericType(value, path);
}

//------------------------------------------------------------------------------
// Read an int from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(int& value, const string path) const {
  readGenericType(value, path);
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(bool& value, const string path) const {
  readGenericType(value, path);
}

//------------------------------------------------------------------------------
// Read a Scalar from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(double& value, const string path) const {
  readGenericType(value, path);
}

//------------------------------------------------------------------------------
// Read a Vector from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(Dim<1>::Vector& value, const string path) const {
  readGenericType(value, path);
}

void
FlatFileIO::read(Dim<2>::Vector& value, const string path) const {
  readGenericType(value, path);
}

void
FlatFileIO::read(Dim<3>::Vector& value, const string path) const {
  readGenericType(value, path);
}

//------------------------------------------------------------------------------
// Read a Tensor from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(Dim<1>::Tensor& value, const string path) const {
  readGenericType(value, path);
}

void
FlatFileIO::read(Dim<2>::Tensor& value, const string path) const {
  readGenericType(value, path);
}

void
FlatFileIO::read(Dim<3>::Tensor& value, const string path) const {
  readGenericType(value, path);
}

//------------------------------------------------------------------------------
// Read a SymTensor from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(Dim<1>::SymTensor& value, const string path) const {
  readGenericType(value, path);
}

void
FlatFileIO::read(Dim<2>::SymTensor& value, const string path) const {
  readGenericType(value, path);
}

void
FlatFileIO::read(Dim<3>::SymTensor& value, const string path) const {
  readGenericType(value, path);
}

//------------------------------------------------------------------------------
// Read a ThirdRankTensor from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(Dim<1>::ThirdRankTensor& value, const string path) const {
  readGenericType(value, path);
}

void
FlatFileIO::read(Dim<2>::ThirdRankTensor& value, const string path) const {
  readGenericType(value, path);
}

void
FlatFileIO::read(Dim<3>::ThirdRankTensor& value, const string path) const {
  readGenericType(value, path);
}

//------------------------------------------------------------------------------
// Read a string from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::read(string& value, const string path) const {
  REQUIRE(readyToRead());

  // Set the current position to the beginning of the file.
//   mFilePtr->seekg(ios::beg);
  beginningOfFile();

  // Find the requested path, and read the value.
  findPathName(path);
  getline(*mFilePtr, value);

  // We pick up the leading space in the string, so get rid of that.
  if (value.size() > 0) value.erase(0,1);
}

//------------------------------------------------------------------------------
// Read a std::vector<int> from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
read(std::vector<int>& value, const string path) const {
  std::string strbuf;
  this->read(strbuf, path);
  const std::vector<char> buf(strbuf.begin(), strbuf.end());
  auto itr = buf.begin();
  value.clear();
  unpackElement(value, itr, buf.end());
  CHECK(itr == buf.end());
}

//------------------------------------------------------------------------------
// Read a std::vector<double> from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
read(std::vector<double>& value, const string path) const {
  std::string strbuf;
  this->read(strbuf, path);
  const std::vector<char> buf(strbuf.begin(), strbuf.end());
  auto itr = buf.begin();
  value.clear();
  unpackElement(value, itr, buf.end());
  CHECK(itr == buf.end());
}

//------------------------------------------------------------------------------
// Read a std::vector<string> from the file.
//------------------------------------------------------------------------------
void
FlatFileIO::
read(std::vector<string>& value, const string path) const {
  std::string strbuf;
  this->read(strbuf, path);
  const std::vector<char> buf(strbuf.begin(), strbuf.end());
  auto itr = buf.begin();
  value.clear();
  unpackElement(value, itr, buf.end());
  CHECK(itr == buf.end());
}

//------------------------------------------------------------------------------
// Access the precision we're using to write to the file.
//------------------------------------------------------------------------------
int
FlatFileIO::
precision() const {
  return mPrecision;
}

void
FlatFileIO::
setPrecision(int precision) {
  REQUIRE(precision > 0);
  mPrecision = precision;
}

//------------------------------------------------------------------------------
// Generic function to write an arbitrary DataType to a FlatFile file.
//------------------------------------------------------------------------------
template<typename DataType>
void
FlatFileIO::
writeGenericType(const DataType& value,
		 const string path) {
  REQUIRE(readyToWrite());
//   mFilePtr->seekp(ios::end);
  *mFilePtr << path << " " << value << endl;
}

//------------------------------------------------------------------------------
// Generic function to read an arbitrary DataType from a FlatFile file.
//------------------------------------------------------------------------------
template<typename DataType>
void
FlatFileIO::
readGenericType(DataType& value,
		const string path) const {


  REQUIRE(readyToRead());

  // Set the current position to the beginning of the file.
//   mFilePtr->seekg(ios::beg);
  beginningOfFile();

  // Find the requested path, and read the value.
  findPathName(path);
  *mFilePtr >> value;

}

//------------------------------------------------------------------------------
// Generic function to write a vector<DataType> to a FlatFile file.
//------------------------------------------------------------------------------
template<typename DataType>
void
FlatFileIO::
writeGenericVector(const vector<DataType>& value,
                   const string path) {
  REQUIRE(readyToWrite());
//   mFilePtr->seekp(ios::end);
  *mFilePtr << path << " " << value.size();
  for (typename vector<DataType>::const_iterator itr = value.begin();
       itr < value.end();
       ++itr) *mFilePtr << " " << *itr;
  *mFilePtr << endl;
}

//------------------------------------------------------------------------------
// Generic function to read a vector<DataType> from a FlatFile file.
//------------------------------------------------------------------------------
template<typename DataType>
void
FlatFileIO::
readGenericVector(vector<DataType>& value,
                  const string path) const {


  REQUIRE(readyToRead());

  // Set the current position to the beginning of the file.
//   mFilePtr->seekg(ios::beg);
  beginningOfFile();

  // Find the requested path, and read the value.
  findPathName(path);

  // Read the size of the stored vector data.
  int size;
  *mFilePtr >> size;
  value.resize(size);

  // Now read in the vector values.
  for (typename vector<DataType>::iterator itr = value.begin();
       itr < value.end();
       ++itr) *mFilePtr >> *itr;

}

//------------------------------------------------------------------------------
// A function to determine if the current file is ready to write to.
//------------------------------------------------------------------------------
bool
FlatFileIO::readyToWrite() const {
  return (mFilePtr != 0 &&
	  (access() == AccessType::Write || access() == AccessType::ReadWrite || access() == AccessType::Create));
}

//------------------------------------------------------------------------------
// A function to determine if the current file is ready to read from.
//------------------------------------------------------------------------------
bool
FlatFileIO::readyToRead() const {
  return (mFilePtr != 0 && 
	  (access() == AccessType::Read || access() == AccessType::ReadWrite || access() == AccessType::Create));
}

//------------------------------------------------------------------------------
// Scan the file until we find the given path name.
//------------------------------------------------------------------------------
void
FlatFileIO::findPathName(const string path) const {
  REQUIRE(readyToRead());


  // Set the pointer to the beginning of the file stream.
//   mFilePtr->seekg(ios::beg);
  beginningOfFile();

  // Scan the file line by line until we find the line containing 
  // the path, or we hit the end.
  string currentPath;
  while (!mFilePtr->eof() && currentPath != path) {
    currentPath = "";
    char thpt = '0';
    size_t i = 0;
    while (!mFilePtr->eof() && thpt != '\n' && i < path.size()) {
      mFilePtr->get(thpt);
      currentPath += thpt;
      ++i;
    }
//     cerr << path << " " << currentPath << " " << (path == currentPath) << endl;

    // If this line doesn't contain the path, skip to the next line.
    if (currentPath != path) {
      while (!mFilePtr->eof() && thpt != '\n') mFilePtr->get(thpt);
    }
  }

  VERIFY2(currentPath == path,
          "FlatFileIO::findPathName ERROR: couldn't find path " << path);
}

//------------------------------------------------------------------------------
// Position the current pointer at the beginning of the file.
//------------------------------------------------------------------------------
void
FlatFileIO::beginningOfFile() const {

//   mFilePtr->seekg(0);
//   ENSURE(mFilePtr->tellg() == 0);

//   close();
//   open(fileName(), access());

  delete mFilePtr;
  mFilePtr = new std::fstream(fileName().c_str(), ios::in);
}

}
