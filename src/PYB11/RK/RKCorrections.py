#-------------------------------------------------------------------------------
# RKCorrections
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class RKCorrections(Physics):
    "Computes RK correction terms"
    
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
               orders = "const std::set<RKOrder>",
               dataBase = "const DataBase<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               volumeType = "const RKVolumeType",
               needHessian = "const bool",
               updateInFinalize = ("const bool", False)):
        "Constructor"
        
    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Evaluate derivatives"
        return "void"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state"
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register derivatives"
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields"
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields"
        return "void"
    
    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "Tasks we do once on problem startup"
        return "void"

    @PYB11virtual
    def initializeProblemStartupDependencies(self,
                                             dataBase = "DataBase<%(Dimension)s>&",
                                             state = "State<%(Dimension)s>&",
                                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Tasks we do once on problem startup"
        return "void"

    @PYB11virtual
    def preStepInitialize(self,
                          dataBase = "const DataBase<%(Dimension)s>&", 
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Optional hook to be called at the beginning of a time step"
        return "void"

    @PYB11virtual
    def initialize(self,
                   time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize the Hydro before we start a derivative evaluation."
        return "bool"
                  
    @PYB11virtual
    def finalize(self,
                 time = "const Scalar",
                 dt = "const Scalar",
                 dataBase = "DataBase<%(Dimension)s>&",
                 state = "State<%(Dimension)s>&",
                 derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize the hydro at the completion of an integration step."
        return "void"
                  
    @PYB11virtual
    def addFacetedBoundary(self,
                           cell = "const FacetedVolume&",
                           holes = ("const std::vector<FacetedVolume>&", "std::vector<FacetedVolume>()")):
        "Add a faceted boundary for Voronoi volume calculation."
        return "void"
        
    @PYB11const
    @PYB11returnpolicy("reference_internal")
    @PYB11keepalive(0,1)
    def WR(self,
           order = "const RKOrder"):
        "Look up the ReproducingKernel for the given order"
        return "const ReproducingKernel<%(Dimension)s>&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    @PYB11keepalive(0,1)
    def corrections(self,
                    order = "const RKOrder"):
        "Look up the corrections for the given order"
        return "const FieldList<%(Dimension)s, RKCoefficients<%(Dimension)s>>&"

    #...........................................................................
    # Properties
    correctionOrders = PYB11property(doc="The set of spatial orders for the reproducing kernel corrections")
    volumeType = PYB11property(doc="Flag for the RK volume weighting definition")
    needHessian = PYB11property(doc="Flag for the RK volume weighting definition")
    volume = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "volume", returnpolicy="reference_internal")

    surfaceArea = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "surfaceArea", returnpolicy="reference_internal")
    normal = PYB11property("const FieldList<%(Dimension)s, Vector>&", "normal", returnpolicy="reference_internal")
    surfacePoint = PYB11property("const FieldList<%(Dimension)s, int>&", "surfacePoint", returnpolicy="reference_internal")
    etaVoidPoints = PYB11property("const FieldList<%(Dimension)s, std::vector<Vector>>&", "etaVoidPoints", returnpolicy="reference_internal")
    cells = PYB11property("const FieldList<%(Dimension)s, FacetedVolume>&", "cells", returnpolicy="reference_internal")
    cellFaceFlags = PYB11property("const FieldList<%(Dimension)s, std::vector<CellFaceFlag>>&", "cellFaceFlags", returnpolicy="reference_internal")
    deltaCentroid = PYB11property("const FieldList<%(Dimension)s, Vector>&", "deltaCentroid", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, RKCorrections, virtual=True)
PYB11inject(RestartMethods, RKCorrections)

