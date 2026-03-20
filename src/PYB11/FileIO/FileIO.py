#-------------------------------------------------------------------------------
# FileIO abstract class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from FileIOAbstractMethods import *
from FileIOTemplateMethods import *
from spheralDimensions import *
dims = spheralDimensions()

@PYB11module("SpheralCompiledModules.SpheralFileIO")
class FileIO:
    "Abstract base class for FileIO objects"

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                filename = "const std::string",
                access = "AccessType"):
        "Construct with a given file name and access"

    #...........................................................................
    # Abstract methods
    @PYB11pure_virtual
    def open(self,
             fileName = "const std::string",
             access = "AccessType"):
        "Open a file for IO"
        return "void"

    @PYB11pure_virtual
    def close(self):
        "Close the current file we're pointing at"
        return "void"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11noconvert
    def write_unsigned_int(self,
                           value = "const unsigned int",
                           pathName = "const std::string"):
        "Write an unsigned int"
        return "void"

    @PYB11virtual
    @PYB11noconvert
    def write_size_t(self,
                     value = "const size_t",
                     pathName = "const std::string"):
        "Write a size_t"
        return "void"

    @PYB11virtual
    @PYB11noconvert
    def write_int(self,
                  value = "const int",
                  pathName = "const std::string"):
        "Write an int"
        return "void"

    @PYB11virtual
    @PYB11noconvert
    def write_bool(self,
                   value = "const bool",
                   pathName = "const std::string"):
        "Write a bool"
        return "void"

    @PYB11virtual
    @PYB11noconvert
    def write_double(self,
                     value = "const double",
                     pathName = "const std::string"):
        "Write a double"
        return "void"

    @PYB11virtual
    @PYB11noconvert
    def write_string(self,
                     value = "const std::string",
                     pathName = "const std::string"):
        "Write a std::string"
        return "void"

    @PYB11virtual
    @PYB11noconvert
    def write_vector_char(self,
                          value = "const std::vector<char>&",
                          pathName = "const std::string"):
        "Write std::vector<char>"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11noconvert
    def read_unsigned_int(self,
                          pathName = "const std::string"):
        "Read an unsigned int"
        return "unsigned int"

    @PYB11virtual
    @PYB11const
    @PYB11noconvert
    def read_size_t(self,
                    pathName = "const std::string"):
        "Read a size_t"
        return "size_t"

    @PYB11virtual
    @PYB11const
    @PYB11noconvert
    def read_int(self,
                 pathName = "const std::string"):
        "Read an int"
        return "int"

    @PYB11virtual
    @PYB11const
    @PYB11noconvert
    def read_bool(self,
                  pathName = "const std::string"):
        "Read a bool"
        return "bool"

    @PYB11virtual
    @PYB11const
    @PYB11noconvert
    def read_double(self,
                    pathName = "const std::string"):
        "Read a double"
        return "double"

    @PYB11virtual
    @PYB11const
    @PYB11noconvert
    def read_string(self,
                    pathName = "const std::string"):
        "Read a std::string"
        return "std::string"

    @PYB11virtual
    @PYB11const
    @PYB11noconvert
    def read_vector_char(self,
                         pathName = "const std::string"):
        "Read a std::vector<char>"
        return "std::vector<char>"

    #...........................................................................
    for ndim in range(1,4):  # These dimensional methods are always supported
        exec('''
@PYB11pycppname("write")
@PYB11noconvert
def writePlane%(ndim)i(self,
                       value = "const GeomPlane<Dim<%(ndim)i>>&",
                       pathName = "const std::string"):
    "Write a Plane%(ndim)id"
    return "void"

@PYB11pycppname("read")
@PYB11const
@PYB11noconvert
def readPlane%(ndim)i(self,
                      value = "GeomPlane<Dim<%(ndim)i>>&",
                      pathName = "const std::string"):
    "Read a Plane%(ndim)id"
    return "void"

@PYB11pycppname("write")
@PYB11noconvert
def writeFacetedVolume%(ndim)i(self,
                               value = "const Dim<%(ndim)i>::FacetedVolume&",
                               pathName = "const std::string"):
    "Write a FacetedVolume%(ndim)i"
    return "void"

@PYB11pycppname("read")
@PYB11const
@PYB11noconvert
def readFacetedVolume%(ndim)i(self,
                              value = "Dim<%(ndim)i>::FacetedVolume&",
                              pathName = "const std::string"):
    "Read a FacetedVolume%(ndim)i"
    return "void"
''' % {"ndim" : ndim})

    @PYB11pycppname("write")
    def write_uniform_random(self,
                             value = "const uniform_random&",
                             pathName = "const std::string"):
        "Write random number generator uniform_random"
        return "void"

    @PYB11pycppname("read")
    @PYB11const
    def read_uniform_random(self,
                            value = "uniform_random&",
                            pathName = "const std::string"):
        "Read random number generator uniform_random"
        return "void"

    @PYB11const
    def splitPathComponents(self, pathName="const std::string"):
        "A helper function to split a string up into substrings delimited by '/'."
        return "std::vector<std::string>"

    @PYB11const
    def groupName(self, pathName="const std::string"):
        "Return the group (directory) component of a path."
        return "std::string"

    @PYB11const
    def variableName(self, pathName="const std::string"):
        "Return the variable component of a path."
        return "std::string"

    def writeObject(self,
                     thing = "py::object",
                     path = "const std::string"):
        "Handle a generic python object through serialization"
        return "void"

    @PYB11returnpolicy("take_ownership")
    @PYB11const
    def readObject(self,
                    path = "const std::string"):
        "Return a generic python object from deserialization."
        return "py::object"

    @PYB11virtual
    def writeBytes(self,
                   stuff = "py::bytes",
                   path = "const std::string"):
        "Specialized method for writing py::bytes -- can be overridden"
        return "void"

    @PYB11virtual
    @PYB11const
    def readBytes(self,
                  path = "const std::string"):
        "Specilized method for reading py::bytes -- can be overridden"
        return "py::bytes"

    @PYB11virtual
    @PYB11pycppname("write")
    def write_py_bytes(self,
                       stuff = "py::bytes&",
                       path = "const std::string"):
        "Override generic write for py::bytes"
        return "void"

    @PYB11virtual
    @PYB11const
    @PYB11pycppname("read")
    def read_py_bytes(self,
                      stuff = "py::bytes&",
                      path = "const std::string"):
        "Override generic read for py::bytes"
        return "void"

    #...........................................................................
    # Properties
    fileName = PYB11property("const std::string&", "fileName", doc="The current file name")
    access = PYB11property("AccessType", "access", doc="The access type of the currently open file")
    fileOpen = PYB11property("bool", "fileOpen", doc="Is the file currently open?")

#-------------------------------------------------------------------------------
# Inject the abstract interface
#-------------------------------------------------------------------------------
PYB11inject(FileIOAbstractMethods, FileIO, virtual=False, pure_virtual=True)
PYB11inject(FileIOTemplateMethods, FileIO)
