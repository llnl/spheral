#-------------------------------------------------------------------------------
# FSISPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralFSISPH")
class SlideSurface:
    "SlideSurface -- some helper methods for to implement SPH slidelines"

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
"""
    
    def pyinit(dataBase = "DataBase<%(Dimension)s>&",
               contactTypes = "const vector<int>"):
        "Slide surface constructor"

    @PYB11const
    def isSlideSurface(nodeListi = "const int",
                       nodeListj = "const int"):
        "returns true if slide surface between nodelisti and j"
        return "bool"

    #...........................................................................
    # Properties
    numNodeLists = PYB11property("int", "numNodeLists", "numNodeLists", doc="number of nodelists.")
    isActive = PYB11property("bool", "isActive", "isActive", doc="switch if we have a slide.")
