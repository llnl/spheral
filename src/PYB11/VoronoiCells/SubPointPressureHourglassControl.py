#-------------------------------------------------------------------------------
# SubPointPressureHourglassControl
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11dynamic_attr
class SubPointPressureHourglassControl(Physics):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using FacetedVolume = typename %(Dimension)s::FacetedVolume;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""
    
    def pyinit(self,
               fHG = "Scalar",
               xfilter = ("Scalar", "0.0")):
        "SubPointPressureHourglassControl constructor"
        return

    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        return "void"

    @PYB11virtual
    @PYB11const
    def requireVoronoiCells(self):
        "Some physics algorithms require the Voronoi cells per point be computed."
        return "bool"

    #...........................................................................
    # Properties
    fHG = PYB11property("Scalar", "fHG", "fHG", doc="The fractional multiplier on the hourglass force")           
    xfilter = PYB11property("Scalar", "xfilter", "xfilter", doc="The fractional multiplier on the hourglass centroidal position filter")           
    DvDt = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DvDt", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, SubPointPressureHourglassControl)
PYB11inject(RestartMethods, SubPointPressureHourglassControl)
