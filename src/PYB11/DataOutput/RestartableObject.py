#-------------------------------------------------------------------------------
# RestartableObject
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11module("SpheralCompiledModules.SpheralDataOutput")
class RestartableObject:
    "The base class for building restartable python objects in Spheral."

    def pyinit(self,
               pyself = "py::handle",
               priority = ("const unsigned", "100")):
        "Construct with the given object and priority."
        return

    @PYB11virtual
    @PYB11const
    def label(self):
        "Define the label for storing this object in a restart file."
        return "std::string"

    @PYB11virtual
    @PYB11const
    def dumpState(self,
                  file = "FileIO&",
                  pathName = "const std::string&"):
        "Write this objects state to the file under the given path."
        return "void"

    @PYB11virtual
    def restoreState(self,
                     file = "const FileIO&",
                     pathName = ("const std::string&", '"pathName"')):
        "Read the state for this object from the given file and path."
        return "void"
