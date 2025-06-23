#-------------------------------------------------------------------------------
# SPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from RestartMethods import *
from Physics import *

@PYB11template("Dimension")
@PYB11module("SpheralDEM")
class DEMBase(Physics):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using RotationType = typename DEMDimension<%(Dimension)s>::AngularVector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""
    
    def pyinit(self,
               dataBase = "const DataBase<%(Dimension)s>&",
               stepsPerCollision = "const Scalar",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "DEMBase constructor"

    #...........................................................................
    # Virtual methods

    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "Tasks we do once on problem startup."
        return "void"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state Hydro expects to use and evolve."
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11virtual
    def preStepInitialize(self,
                          dataBase = "const DataBase<%(Dimension)s>&", 
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Optional hook to be called at the beginning of a time step."
        return "void"

    @PYB11virtual
    def initialize(self,
                   time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize the DEM before we start a derivative evaluation."
        return "bool"                

    @PYB11virtual
    @PYB11const
    def finalizeDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize the derivatives."
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    @PYB11virtual
    def resizePairFieldLists(self):
        "resize all pair fieldlists consistent w/ neighborIndices"
        return "void"

    @PYB11virtual
    @PYB11const
    def resizeStatePairFieldLists(self,
                                  state = "State<%(Dimension)s>&"):
        "resize pair fieldlists belonging to the state consistent w/ neighborIndices"
        return "void"

    @PYB11virtual
    @PYB11const
    def resizeDerivativePairFieldLists(self,
                                       derivs = "StateDerivatives<%(Dimension)s>&"):
        "resize pair fieldlists belonging to the derivatives consistent w/ neighborIndices"
        return "void"

    @PYB11virtual
    def removeInactiveContactsFromPairFieldLists(self):
        "remove old contacts from all pair fieldlists "
        return "void"

    @PYB11virtual
    @PYB11const
    def removeInactiveContactsFromStatePairFieldLists(self,
                                  state = "State<%(Dimension)s>&"):
        "remove old contacts from pair fieldlists belonging to the state consistent w/ neighborIndices"
        return "void"

    @PYB11virtual
    @PYB11const
    def removeInactiveContactsFromDerivativePairFieldLists(self,
                                       derivs = "StateDerivatives<%(Dimension)s>&"):
        "remove old contacts from pair fieldlists belonging to the derivatives consistent w/ neighborIndices"
        return "void"

    def initializeOverlap(self,
                          dataBase = "const DataBase<%(Dimension)s>&",
                          startCompositeParticleIndex = "const int"):
        "set the equilibrium overlap pairwise fieldlist for comp. particle id's > specified value"
        return "void"

    def updateContactMap(self,
                         dataBase = "const DataBase<%(Dimension)s>&"):
        "update DEM contact/neighbor tracker"
        return "void"
    
    def identifyInactiveContacts(self,
                                 dataBase = "const DataBase<%(Dimension)s>&"):
        "initializes the isActive pairfieldlist to flag outdated contacts"
        return "void"

    def updatePairwiseFieldLists(self,
                                 purgeInactiveContacts = "const bool"):
        "update DEM contact map and pair fieldlists consistent w/ the connectivityMap"
        return "void"

    def appendSolidBoundary(self,
                            boundary = "SolidBoundaryBase<%(Dimension)s>&"):
        "add a solid boundary to the end of the list"
        return "void"

    def clearSolidBoundaries(self):
        "remove all solid boundaries from the dem package"
        return "void"

    @PYB11const
    def haveSolidBoundary(self,
                          boundary = "const SolidBoundaryBase<%(Dimension)s>&"):
        "is this boundary being tracked?"
        return "bool"
    
    def removeSolidBoundary(self,
                            boundary = "const SolidBoundaryBase<%(Dimension)s>&"):
        "remove the specified solid boundary"
        return "void"

    #...........................................................................
    # Properties
    xmin = PYB11property("const Vector&", "xmin", "xmin",
                         returnpolicy="reference_internal",
                         doc="Optional minimum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    xmax = PYB11property("const Vector&", "xmax", "xmax",
                         returnpolicy="reference_internal",
                         doc="Optional maximum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    
    stepsPerCollision = PYB11property("Scalar", "stepsPerCollision", "stepsPerCollision", doc="number of steps resolving the collision time scale")
    
    cycle = PYB11property("int", "cycle", "cycle", doc="tracks what cycle we are on.")
    contactRemovalFrequency = PYB11property("int", "contactRemovalFrequency", "contactRemovalFrequency", doc="every n-cycles we prune our contact list to only the active contacts. This is the frequency.")

    timeStepMask =  PYB11property("const FieldList<%(Dimension)s, int>&", "timeStepMask", returnpolicy="reference_internal")
    DxDt =          PYB11property("const FieldList<%(Dimension)s, Vector>&","DxDt", returnpolicy="reference_internal")
    DvDt =          PYB11property("const FieldList<%(Dimension)s, Vector>&", "DvDt", returnpolicy="reference_internal")
    DomegaDt =      PYB11property("const FieldList<%(Dimension)s, RotationType>&","DomegaDt", returnpolicy="reference_internal")
    omega =         PYB11property("const FieldList<%(Dimension)s, RotationType>&","omega", returnpolicy="reference_internal")

    equilibriumOverlap = PYB11property("const FieldList<%(Dimension)s, vector<Scalar>>&","equilibriumOverlap", returnpolicy="reference_internal")
    neighborIndices = PYB11property("const FieldList<%(Dimension)s, vector<int>>&","neighborIndices", returnpolicy="reference_internal")
    shearDisplacement = PYB11property("const FieldList<%(Dimension)s, vector<Vector>>&","shearDisplacement", returnpolicy="reference_internal")
    DDtShearDisplacement = PYB11property("const FieldList<%(Dimension)s, vector<Vector>>&","DDtShearDisplacement", returnpolicy="reference_internal")
    isActiveContact = PYB11property("const FieldList<%(Dimension)s, vector<int>>&","isActiveContact", returnpolicy="reference_internal")
    
    newSolidBoundaryIndex = PYB11property("int", "newSolidBoundaryIndex", doc="index of the most recent solid bc added to the package")
    numSolidBoundaries = PYB11property("unsigned int", "numSolidBoundaries", doc="number of solid boundaries")
    numContacts = PYB11property("unsigned int", "numContacts", doc="Total number of contacts")
    numParticleParticleContacts = PYB11property("unsigned int", "numParticleParticleContacts", doc="Number of interactions with other dem particles")
    numParticleBoundaryContacts = PYB11property("unsigned int", "numParticleBoundaryContacts", doc="Number interactions with solid boundaries")
    contactStorageIndices = PYB11property("const std::vector<ContactIndex>&", "contactStorageIndices")
    solidBoundaryConditions = PYB11property("const std::vector<SolidBoundaryBase<%(Dimension)s>*>&", "solidBoundaryConditions", doc="The set of NodeLists in the DataBase")
#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, DEMBase)

