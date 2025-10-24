include_directories(.)
set(FileIO_sources
  FileIO.cc
  FlatFileIO.cc
  SidreFileIO.cc
  SiloFileIO.cc
  vectorstringUtilities.cc)

set(FileIO_headers
  FileIO.hh
  FileIOInline.hh
  FlatFileIO.hh
  HDF5IO.hh
  HDF5Traits.hh
  HDF5Types.hh
  SidreFileIO.hh
  SiloFileIO.hh
  vectorstringUtilities.hh)

if (SPHERAL_ENABLE_PYTHON)
  list(APPEND FileIO_sources PyFileIO.cc)
  list(APPEND FileIO_headers PyFileIO.hh)
endif()

spheral_install_python_files(
  GzipFileIO.py
  PlyFileIO.py)

spheral_add_obj_library(FileIO SPHERAL_OBJ_LIBS)

