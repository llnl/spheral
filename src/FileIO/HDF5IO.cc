//---------------------------------Spheral++----------------------------------//
// HDF5IO -- Provide the interface to HDF5 file objects.
//
// Created by JMO, Wed Jul 12 17:04:37 PDT 2000
//----------------------------------------------------------------------------//

#include "HDF5IO.hh"
#include "Utilities/DBC.hh"
#include "Field/Field.hh"
#include "Utilities/Logger.hh"

#include <sstream>

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
HDF5IO<Dimension>::HDF5IO():
  FileIO<Dimension>(),
  mFilePtr(0) {

  DEBUG_LOG << "HDF5IO::HDF5IO()";

  initializeAccessMap();

  // Turn off error printing from HDF5 by default.
  FileIException error;
  error.dontPrint();

}

//------------------------------------------------------------------------------
// Construct and open the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
HDF5IO<Dimension>::HDF5IO(const string& fileName, AccessType access):
  FileIO<Dimension>(fileName, access),
  mFilePtr(0) {
  DEBUG_LOG << "HDF5IO::HDF5IO(string, access)";
  initializeAccessMap();

  // Turn off error printing from HDF5 by default.
  FileIException error;
  error.dontPrint();

  open(fileName, access);
  CHECK(mFileOpen && mFilePtr != 0);
}

// Do the same thing with an int in the place of the enum, to help Pyffle.
template<typename Dimension>
HDF5IO<Dimension>::HDF5IO(const string& fileName, int access):
  FileIO<Dimension>(fileName, AccessType(access)),
  mFilePtr(0) {
  DEBUG_LOG << "HDF5IO::HDF5IO(string, int)";
  initializeAccessMap();

  // Turn off error printing from HDF5 by default.
  FileIException error;
  error.dontPrint();

  open(fileName, AccessType(access));
  CHECK(mFileOpen && mFilePtr != 0);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
HDF5IO<Dimension>::~HDF5IO() {
  DEBUG_LOG << "HDF5IO::~HDF5IO";
  close();
}

//------------------------------------------------------------------------------
// Open an HDF5 file with the specified access.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::open(const string& fileName, AccessType access) {

  DEBUG_LOG << "HDF5IO::open(" << fileName << ", " << access << ")";

  // If we currently have a file open and attached to this object, close it!
  close();
  CHECK(mFilePtr == 0 && mFileOpen == false);

  // Open the new file.
  try {
    mFilePtr = new H5File(fileName, mHDF5AccessTypes[access]);
    mFileOpen = true;
  } catch(FileIException fileError) {
    fileError.printError();
  }
  CHECK(mFilePtr != 0 && mFileOpen == true);
}

//------------------------------------------------------------------------------
// Close the current file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::close() {

  DEBUG_LOG << "HDF5IO::close()";

  if (mFilePtr != 0) {
    delete mFilePtr;
    mFilePtr = 0;
  }
  mFileOpen = false;
}

//------------------------------------------------------------------------------
// Write an int to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::write(const int& value, const string& pathName) {
  writeGenericType(value, pathName, PredType::NATIVE_INT);
}

//------------------------------------------------------------------------------
// Write a bool to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::write(const bool& value, const string& pathName) {
  writeGenericType(value, pathName, PredType::NATIVE_HBOOL);
}

//------------------------------------------------------------------------------
// Write a Scalar to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
write(const typename Dimension::Scalar& value, const string& pathName) {
  writeGenericType(value, pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Write a Vector to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
write(const typename Dimension::Vector& value, const string& pathName) {
  writeGenericContainer(value.begin(), value.end(), pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Write a Tensor to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
write(const typename Dimension::Tensor& value, const string& pathName) {
  writeGenericContainer(value.begin(), value.end(), pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Write a SymTensor to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
write(const typename Dimension::SymTensor& value, const string& pathName) {
  writeGenericContainer(value.begin(), value.end(), pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Write a string to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
write(const string& value, const string& pathName) {
  writeGenericContainer(value.begin(), value.end(), pathName, PredType::NATIVE_CHAR);
}

//------------------------------------------------------------------------------
// Read an int from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::read(int& value, const string& pathName) const {
  readGenericType(value, pathName, PredType::NATIVE_INT);
}

//------------------------------------------------------------------------------
// Read a bool from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::read(bool& value, const string& pathName) const {
  readGenericType(value, pathName, PredType::NATIVE_HBOOL);
}

//------------------------------------------------------------------------------
// Read a Scalar from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::read(typename Dimension::Scalar& value, const string& pathName) const {
  readGenericType(value, pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Read a Vector from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::read(typename Dimension::Vector& value, const string& pathName) const {
  readGenericContainer(value.begin(), value.end(), pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Read a Tensor from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::read(typename Dimension::Tensor& value, const string& pathName) const {
  readGenericContainer(value.begin(), value.end(), pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Read a SymTensor from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::read(typename Dimension::SymTensor& value, const string& pathName) const {
  readGenericContainer(value.begin(), value.end(), pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Read a string from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::read(string& value, const string& pathName) const {
  readGenericContainer(value.begin(), value.end(), pathName, PredType::NATIVE_CHAR);
}

//------------------------------------------------------------------------------
// Write a Scalar Field to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
write(const Field<Dimension, typename Dimension::Scalar>& value, const string& pathName) {
  writeGenericContainer(value.internalBegin(), value.internalEnd(), pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Write a Vector Field to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
write(const Field<Dimension, typename Dimension::Vector>& value, const string& pathName) {
  writeGenericField(value, pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Write a Tensor Field to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
write(const Field<Dimension, typename Dimension::Tensor>& value, const string& pathName) {
  writeGenericField(value, pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Write a SymTensor Field to the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
write(const Field<Dimension, typename Dimension::SymTensor>& value, const string& pathName) {
  writeGenericField(value, pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Read a Scalar Field from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
read(Field<Dimension, typename Dimension::Scalar>& value, const string& pathName) const {
  readGenericContainer(value.internalBegin(), value.internalEnd(), pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Read a Vector Field from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
read(Field<Dimension, typename Dimension::Vector>& value, const string& pathName) const {
  readGenericField(value, pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Read a Tensor Field from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
read(Field<Dimension, typename Dimension::Tensor>& value, const string& pathName) const {
  readGenericField(value, pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Read a SymTensor Field from the file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::
read(Field<Dimension, typename Dimension::SymTensor>& value, const string& pathName) const {
  readGenericField(value, pathName, PredType::NATIVE_DOUBLE);
}

//------------------------------------------------------------------------------
// Generic function to write an arbitrary DataType to an HDF5 file.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType, typename H5DataType>
void
HDF5IO<Dimension>::
writeGenericType(const DataType& value,
		 const string& pathName,
		 const H5DataType& h5Type) {

  CHECK(readyToWrite());

  // Get the group portion of the path name.
  string grpName = groupName(pathName);
  CHECK(grpName != "");

  // Get the variable name we are using for this variable.
  string varName = variableName(pathName);
  CHECK(varName != "");

  // Catch exceptions thrown by HDF5.
  try {

    // Create a DataSpace to hold a single element.
    hsize_t dim[1] = {1};
    DataSpace dataSpace(1, dim);

    // Get the group representing this path name.
    Group* groupPtr = getH5Group(grpName);
    CHECK(groupPtr != 0);

    // Create the DataSet.
    DataSet dataSet = groupPtr->createDataSet(varName, h5Type, dataSpace);

    // Write the data to the file.
    dataSet.write(&value, h5Type);
    //    dataSet.flush(H5F_SCOPE_LOCAL);

    // Release the group.
    //    groupPtr->flush(H5F_SCOPE_LOCAL);
    delete groupPtr;
  }

  // catch failure caused by the H5File operations
  catch(FileIException error) {
    error.printError();
  }

  // catch failure caused by the DataSet operations
  catch(DataSetIException error) {
    error.printError();
  }

  // catch failure caused by the DataSpace operations
  catch(DataSpaceIException error) {
    error.printError();
  }

}

//------------------------------------------------------------------------------
// Generic function to write a container of arbitrary data to an HDF5 file.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename IteratorType, typename H5DataType>
void
HDF5IO<Dimension>::
writeGenericContainer(const IteratorType& begin, 
                      const IteratorType& end,
		      const string& pathName,
                      const H5DataType& h5Type) {

  CHECK(readyToWrite());

  // Get the group portion of the path name.
  string grpName = groupName(pathName);
  CHECK(grpName != "");

  // Get the variable name we are using for this variable.
  string varName = variableName(pathName);
  CHECK(varName != "");

  // Catch exceptions thrown by HDF5.
  try {

    // Create a DataSpace to hold a single element.
    hsize_t dim[1] = {distance(begin, end)};
    DataSpace dataSpace(1, dim);

    // Get the group representing this path name.
    Group* groupPtr = getH5Group(grpName);
    CHECK(groupPtr != 0);

    // Create the DataSet.
    DataSet dataSet = groupPtr->createDataSet(varName, h5Type, dataSpace);

    // Write the data to the file.
    dataSet.write(&(*begin), h5Type);
    //    dataSet.flush(H5F_SCOPE_LOCAL);

    // Release the group.
    //    groupPtr->flush(H5F_SCOPE_LOCAL);
    delete groupPtr;
  }

  // catch failure caused by the H5File operations
  catch(FileIException error) {
    error.printError();
  }

  // catch failure caused by the DataSet operations
  catch(DataSetIException error) {
    error.printError();
  }

  // catch failure caused by the DataSpace operations
  catch(DataSpaceIException error) {
    error.printError();
  }

}

//------------------------------------------------------------------------------
// Generic function to write a Field of arbitrary data to an HDF5 file.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType, typename H5DataType>
void
HDF5IO<Dimension>::
writeGenericField(const Field<Dimension, DataType>& field,
		  const string& pathName,
                  const H5DataType& h5Type) {

  CHECK(readyToWrite());

  // Loop over each element of the Field.
  for (Field<Dimension, DataType>::const_iterator elementItr = field.internalBegin();
       elementItr < field.internalEnd();
       ++elementItr) {

    // We will write this element under it's own unique name in the Fields path.
    int elementID = distance(field.internalBegin(), elementItr);
    ostrstream varPath;
    varPath << pathName << "/element" << elementID << '\0';
    string varPathStr = varPath.str();

    // Now write this element to the file.
    writeGenericContainer(elementItr->begin(), elementItr->end(),
			  varPathStr, h5Type);
                          
  }
}

//------------------------------------------------------------------------------
// Generic function to read an atomic type from an HDF5 file.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType, typename H5DataType>
void
HDF5IO<Dimension>::
readGenericType(DataType& value,
		const string& pathName,
		const H5DataType& h5Type) const {

  CHECK(readyToRead());

  // Get the group portion of the path name.
  string grpName = groupName(pathName);
  CHECK(grpName != "");

  // Get the variable name we are using for this variable.
  string varName = variableName(pathName);
  CHECK(varName != "");

  // Catch exceptions thrown by HDF5.
  try {
    // Get the group representing this path name.
    Group group = mFilePtr->openGroup(grpName);

    // Get the DataSet.
    DataSet dataSet = group.openDataSet(varName);

    // Read the data value.
    dataSet.read(&value, h5Type);
  }

  // catch failure caused by the H5File operations
  catch(FileIException error) {
    error.printError();
  }

  // catch failure caused by the DataSet operations
  catch(DataSetIException error) {
    error.printError();
  }

  // catch failure caused by the DataSpace operations
  catch(DataSpaceIException error) {
    error.printError();
  }
}

//------------------------------------------------------------------------------
// Generic function to read a container of data from an HDF5 file.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename IteratorType, typename H5DataType>
void
HDF5IO<Dimension>::
readGenericContainer(const IteratorType& begin,
		     const IteratorType& end,
		     const string& pathName,
		     const H5DataType& h5Type) const {

  CHECK(readyToRead());

  // Get the group portion of the path name.
  string grpName = groupName(pathName);
  CHECK(grpName != "");

  // Get the variable name we are using for this variable.
  string varName = variableName(pathName);
  CHECK(varName != "");

  // Catch exceptions thrown by HDF5.
  try {
    // Get the group representing this path name.
    Group group = mFilePtr->openGroup(grpName);

    // Get the DataSet.
    DataSet dataSet = group.openDataSet(varName);

    // Get the size of the DataSet.
    DataSpace dataSpace = dataSet.getSpace();
    hsize_t dims[1];
    int rank = dataSpace.getSimpleExtentDims(dims);

    // Verify that the container is the correct size!
    if (distance(begin, end) != (int) dims[0]) {
      std::cerr << "HDF5IO::readGenericContainer ERROR: input container wrong size."
                << std::endl
                << "  Container size: " << distance(begin, end) << std::endl
                << "    DataSet size: " << dims[1] << std::endl;
      return;
    }

    // Read the data value.
    dataSet.read(&(*begin), h5Type);
  }

  // catch failure caused by the H5File operations
  catch(FileIException error) {
    error.printError();
  }

  // catch failure caused by the DataSet operations
  catch(DataSetIException error) {
    error.printError();
  }

  // catch failure caused by the DataSpace operations
  catch(DataSpaceIException error) {
    error.printError();
  }
}

//------------------------------------------------------------------------------
// Generic function to read a Field of arbitrary data from an HDF5 file.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType, typename H5DataType>
void
HDF5IO<Dimension>::
readGenericField(Field<Dimension, DataType>& field,
		 const string& pathName,
		 const H5DataType& h5Type) const {

  CHECK(readyToRead());

  // Loop over each element of the Field.
  for (Field<Dimension, DataType>::iterator elementItr = field.internalBegin();
       elementItr < field.internalEnd();
       ++elementItr) {

    // We will write this element under it's own unique name in the Fields path.
    int elementID = distance(field.internalBegin(), elementItr);
    ostrstream varPath;
    varPath << pathName << "/element" << elementID << '\0';
    string varPathStr = varPath.str();

    // Now write this element to the file.
    readGenericContainer(elementItr->begin(), elementItr->end(),
			 varPathStr, h5Type);
                          
  }
}

//------------------------------------------------------------------------------
// Get a pointer to the group at the end of the given path, creating the path
// if necessary.
//------------------------------------------------------------------------------
template<typename Dimension>
Group*
HDF5IO<Dimension>::getH5Group(const string& pathName) {
  CHECK(readyToWrite());

  // Break up the path name into the individual group components.
  vector<string> groupNames = splitPathComponents(pathName);

  // Loop over the group names, and create new groups for each level if necessary.
  string cumulativePath = "";
  for (vector<string>::const_iterator nameItr = groupNames.begin();
       nameItr < groupNames.end();
       ++nameItr) {

    cumulativePath += *nameItr;

    // Check if this much of the group path exists.
    try {
      // Try to read the group from the file if it already exists.
      Group group = mFilePtr->openGroup(cumulativePath);
    }

    // Catch the case that this group does not yet exist, and create it.
    catch(FileIException error) {
      DEBUG_LOG << "Caught file exception, creating group " << cumulativePath;
      Group group = mFilePtr->createGroup(cumulativePath);
    }
  }

  // Finally, open up the ultimate group name and return a pointer.
  return new Group(mFilePtr->openGroup(pathName));
}

//------------------------------------------------------------------------------
// A function to determine if the current file is ready to write to.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
HDF5IO<Dimension>::readyToWrite() const {
  return (mFilePtr != 0 && mFileOpen &&
	  (access() == Write || access() == ReadWrite) || access() == Create);
}

//------------------------------------------------------------------------------
// A function to determine if the current file is ready to read from.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
HDF5IO<Dimension>::readyToRead() const {
  return (mFilePtr != 0 && mFileOpen);
}

//------------------------------------------------------------------------------
// A method to initialize the map of FileIO access types and HDF5 access types.
//------------------------------------------------------------------------------
template<typename Dimension>
void
HDF5IO<Dimension>::initializeAccessMap() {
  DEBUG_LOG << "HDF5IO::initializeAccessMap()";
  mHDF5AccessTypes[Spheral::FileIO::Create] = H5F_ACC_TRUNC;
  mHDF5AccessTypes[Spheral::FileIO::Read] = H5F_ACC_RDONLY;
  mHDF5AccessTypes[Spheral::FileIO::Write] = H5F_ACC_RDWR;
  mHDF5AccessTypes[Spheral::FileIO::ReadWrite] = H5F_ACC_RDWR;
}

}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
template class Spheral::FileIO::HDF5IO< Dim<1> >;
template class Spheral::FileIO::HDF5IO< Dim<2> >;
template class Spheral::FileIO::HDF5IO< Dim<3> >;
