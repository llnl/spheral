#-------------------------------------------------------------------------------
# VoronoiCells
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11dynamic_attr
class VoronoiCells(Physics):

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
               kernelExtent = "const Scalar",
               facetedBoundaries = ("const std::vector<FacetedVolume>&", "std::vector<FacetedVolume>()"),
               facetedHoles = ("const std::vector<std::vector<FacetedVolume>>&", "std::vector<std::vector<FacetedVolume>>()")):
        "VoronoiCells constructor (with C++ types)"
        return

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        """An optional hook to initialize once when the problem is starting up.
This is called after the materials and NodeLists are created. This method
should set the sizes of all arrays owned by the physics package and initialize
independent variables.
It is assumed after this method has been called it is safe to call
Physics::registerState to create full populated State objects."""
        return "void"

    @PYB11virtual
    def initializeProblemStartupDependencies(self,
                                             dataBase = "DataBase<%(Dimension)s>&",
                                             state = "State<%(Dimension)s>&",
                                             derivs = "StateDerivatives<%(Dimension)s>&"):
        """A second optional method to be called on startup, after Physics::initializeProblemStartup has
been called.
One use for this hook is to fill in dependendent state using the State object, such as
temperature or pressure."""
        return "void"

    @PYB11virtual 
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Increment the derivatives."
        return "void"

    @PYB11virtual
    def preStepInitialize(dataBase = "const DataBase<%(Dimension)s>&",
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize at the beginning of a step."
        return "void"

    @PYB11virtual 
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual 
    def registerState(dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state Hydro expects to use and evolve."
        return "void"

    @PYB11virtual
    def registerDerivatives(dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    @PYB11virtual
    def postStateUpdate(time = "const Scalar",
                        dt = "const Scalar",
                        dataBase = "const DataBase<%(Dimension)s>&",
                        state = "State<%(Dimension)s>&",
                        derivs = "StateDerivatives<%(Dimension)s>&"):
        "Provide a hook to be called after the state has been updated and boundary conditions have been enforced."
        return "bool"

    @PYB11virtual
    def addFacetedBoundary(bound = "const FacetedVolume&",
                           holes = "const std::vector<FacetedVolume>&"):
        "Add a faceted boundary (optionally with holes)"
        return "void"

    @PYB11virtual
    @PYB11const
    def requireConnectivity(self):
        "Returns True, we do need connectivity"
        return "bool"

    #...........................................................................
    # Properties
    kernelExtent =                 PYB11property("Scalar",                                                     "kernelExtent", doc="The kernel extent in eta")           
    volume =                       PYB11property("const FieldList<%(Dimension)s, Scalar>&",                    "volume",            returnpolicy="reference_internal")
    weight =                       PYB11property("const FieldList<%(Dimension)s, Scalar>&",                    "weight",            returnpolicy="reference_internal")
    surfacePoint =                 PYB11property("const FieldList<%(Dimension)s, int>&",                       "surfacePoint",      returnpolicy="reference_internal")
    etaVoidPoints =                PYB11property("const FieldList<%(Dimension)s, std::vector<Vector>>&",       "etaVoidPoints",     returnpolicy="reference_internal")
    cells =                        PYB11property("const FieldList<%(Dimension)s, FacetedVolume>&",             "cells",             returnpolicy="reference_internal")
    cellFaceFlags =                PYB11property("const FieldList<%(Dimension)s, std::vector<CellFaceFlag>>&", "cellFaceFlags",     returnpolicy="reference_internal")
    deltaCentroid =                PYB11property("const FieldList<%(Dimension)s, Vector>&",                    "deltaCentroid",     returnpolicy="reference_internal")
    facetedBoundaries =            PYB11property("const std::vector<FacetedVolume>&",                          "facetedBoundaries", returnpolicy="reference_internal")
    facetedHoles =                 PYB11property("const std::vector<std::vector<FacetedVolume>>&",             "facetedHoles",      returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, VoronoiCells)
