#-------------------------------------------------------------------------------
# PyFileIO
#-------------------------------------------------------------------------------
from PYB11Generator import *
from FileIO import *
from FileIOAbstractMethods import *
from FileIOTemplateMethods import *
from spheralDimensions import *
dims = spheralDimensions()

@PYB11module("SpheralCompiledModules.SpheralFileIO")
class PyFileIO(FileIO):
    "PyFileIO -- A python friendly version of the FileIO interface, for use creating python FileIO objects."

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                filename = "const std::string",
                access = "AccessType"):
        "Construct a PyFileIO object with a file name and access"

    #...........................................................................
    # Abstract interface
    @PYB11virtual
    def write_unsigned_int(self,
                           value = "const unsigned int",
                           pathName = "const std::string"):
        "Write an unsigned int"
        return "void"

    @PYB11virtual
    def write_size_t(self,
                     value = "const size_t",
                     pathName = "const std::string"):
        "Write an size_t"
        return "void"

    @PYB11virtual
    def write_int(self,
                  value = "const int",
                  pathName = "const std::string"):
        "Write an int"
        return "void"

    @PYB11virtual
    def write_bool(self,
                   value = "const bool",
                   pathName = "const std::string"):
        "Write a bool"
        return "void"

    @PYB11virtual
    def write_double(self,
                     value = "const double",
                     pathName = "const std::string"):
        "Write a double"
        return "void"

    @PYB11virtual
    def write_string(self,
                     value = "const std::string",
                     pathName = "const std::string"):
        "Write a std::string"
        return "void"

    @PYB11virtual
    def write_vector_char(self,
                          value = "const std::vector<char>&",
                          pathName = "const std::string"):
        "Write a std::vector<char>"
        return "void"

    @PYB11virtual
    def write_vector_int(self,
                         value = "const std::vector<int>&",
                         pathName = "const std::string"):
        "Write a std::vector<int>"
        return "void"

    @PYB11virtual
    def write_vector_double(self,
                            value = "const std::vector<double>&",
                            pathName = "const std::string"):
        "Write a std::vector<double>"
        return "void"

    @PYB11virtual
    def write_vector_string(self,
                            value = "const std::vector<std::string>&",
                            pathName = "const std::string"):
        "Write a std::vector<string>"
        return "void"

    @PYB11virtual
    def write_Vector1d(self,
                       value = "const Dim<1>::Vector&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11virtual
    def write_Tensor1d(self,
                       value = "const Dim<1>::Tensor&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11virtual
    def write_SymTensor1d(self,
                          value = "const Dim<1>::SymTensor&",
                          pathName = "const std::string"):
        return "void"
  
    @PYB11virtual
    def write_ThirdRankTensor1d(self,
                                value = "const Dim<1>::ThirdRankTensor&",
                                pathName = "const std::string"):
        return "void"

    @PYB11virtual
    def write_Vector2d(self,
                       value = "const Dim<2>::Vector&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11virtual
    def write_Tensor2d(self,
                       value = "const Dim<2>::Tensor&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11virtual
    def write_SymTensor2d(self,
                          value = "const Dim<2>::SymTensor&",
                          pathName = "const std::string"):
        return "void"
  
    @PYB11virtual
    def write_ThirdRankTensor2d(self,
                                value = "const Dim<2>::ThirdRankTensor&",
                                pathName = "const std::string"):
        return "void"

    @PYB11virtual
    def write_Vector3d(self,
                       value = "const Dim<3>::Vector&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11virtual
    def write_Tensor3d(self,
                       value = "const Dim<3>::Tensor&",
                       pathName = "const std::string"):
        return "void"
  
    @PYB11virtual
    def write_SymTensor3d(self,
                          value = "const Dim<3>::SymTensor&",
                          pathName = "const std::string"):
        return "void"
  
    @PYB11virtual
    def write_ThirdRankTensor3d(self,
                                value = "const Dim<3>::ThirdRankTensor&",
                                pathName = "const std::string"):
        return "void"
  
    @PYB11virtual
    @PYB11const
    def read_unsigned_int(self,
                          pathName = "const std::string"):
        "Read an unsigned int"
        return "unsigned int"

    @PYB11virtual
    @PYB11const
    def read_size_t(self,
                    pathName = "const std::string"):
        "Read an size_t"
        return "size_t"

    @PYB11virtual
    @PYB11const
    def read_int(self,
                 pathName = "const std::string"):
        "Read an int"
        return "int"

    @PYB11virtual
    @PYB11const
    def read_bool(self,
                  pathName = "const std::string"):
        "Read a bool"
        return "bool"

    @PYB11virtual
    @PYB11const
    def read_double(self,
                    pathName = "const std::string"):
        "Read a double"
        return "double"

    @PYB11virtual
    @PYB11const
    def read_string(self,
                    pathName = "const std::string"):
        "Read a std::string"
        return "std::string"

    @PYB11virtual
    @PYB11const
    def read_vector_char(self,
                         pathName = "const std::string"):
        "Read a std::vector<char>"
        return "std::vector<char>"

    @PYB11virtual
    @PYB11const
    def read_vector_int(self,
                        value = "std::vector<int>&",
                        pathName = "const std::string"):
        "Read a std::vector<int>"
        return "void"

    @PYB11virtual
    @PYB11const
    def read_vector_double(self,
                           value = "std::vector<double>&",
                           pathName = "const std::string"):
        "Read a std::vector<double>&"
        return "void"

    @PYB11virtual
    @PYB11const
    def read_vector_string(self,
                           value = "std::vector<std::string>&",
                           pathName = "const std::string"):
        "Read a std::vector<string>"
        return "void"

    @PYB11virtual
    @PYB11const
    def read_Vector1d(self,
                      pathName = "const std::string"):
        return "Dim<1>::Vector"
  
    @PYB11virtual
    @PYB11const
    def read_Tensor1d(self,
                      pathName = "const std::string"):
        return "Dim<1>::Tensor"
  
    @PYB11virtual
    @PYB11const
    def read_SymTensor1d(self,
                         pathName = "const std::string"):
        return "Dim<1>::SymTensor"
  
    @PYB11virtual
    @PYB11const
    def read_ThirdRankTensor1d(self,
                               pathName = "const std::string"):
        return "Dim<1>::ThirdRankTensor"
  
    @PYB11virtual
    @PYB11const
    def read_Vector2d(self,
                      pathName = "const std::string"):
        return "Dim<2>::Vector"
  
    @PYB11virtual
    @PYB11const
    def read_Tensor2d(self,
                      pathName = "const std::string"):
        return "Dim<2>::Tensor"
  
    @PYB11virtual
    @PYB11const
    def read_SymTensor2d(self,
                         pathName = "const std::string"):
        return "Dim<2>::SymTensor"
  
    @PYB11virtual
    @PYB11const
    def read_ThirdRankTensor2d(self,
                               pathName = "const std::string"):
        return "Dim<2>::ThirdRankTensor"
  
    @PYB11virtual
    @PYB11const
    def read_Vector3d(self,
                      pathName = "const std::string"):
        return "Dim<3>::Vector"
  
    @PYB11virtual
    @PYB11const
    def read_Tensor3d(self,
                      pathName = "const std::string"):
        return "Dim<3>::Tensor"
  
    @PYB11virtual
    @PYB11const
    def read_SymTensor3d(self,
                         pathName = "const std::string"):
        return "Dim<3>::SymTensor"
  
    @PYB11virtual
    @PYB11const
    def read_ThirdRankTensor3d(self,
                               pathName = "const std::string"):
        return "Dim<3>::ThirdRankTensor"
  
#-------------------------------------------------------------------------------
# Override the required virtual interface
#-------------------------------------------------------------------------------
PYB11inject(FileIOAbstractMethods, PyFileIO, virtual=True, pure_virtual=False)
#PYB11inject(FileIOTemplateMethods, PyFileIO)
