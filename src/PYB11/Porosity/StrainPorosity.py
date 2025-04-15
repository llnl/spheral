#-------------------------------------------------------------------------------
# StrainPorosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from PorosityModel import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralPorosity")
class StrainPorosity(PorosityModel):
    """An implementation of strain-alpha porosity model described in two papers:

Wunnemann, Collins, & Melosh, Icarus, 180, 514-527 (2006)
"A strain-based porosity model for use in hydrocode simulations of impacts
 and implications for transient crater growth in porous targets"

Collins, G. S., Melosh, H. J., & Wunnemann, K. (2011).
"Improvements to the epsilon-alpha porous compaction model for simulating impacts into
high-porosity solar system objects.
International Journal of Impact Engineering, 38(6), 434-439.
http://doi.org/10.1016/j.ijimpeng.2010.10.013
"""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using ScalarField = Field<%(Dimension)s, Scalar>;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "const SolidNodeList<%(Dimension)s>&",
               phi0 = "const double",
               epsE = "const double",
               epsX = "const double",
               kappa = "const double",
               gammaS0 = "const double",
               cS0 = "const double",
               c0 = "const double",
               rhoS0 = "const double",
               jutziStateUpdate = "const bool"):
        """Constructor parameters:
        nodeList:        The NodeList we're going apply to
        phi0:            Initial porosity (single value)
        epsE:            Elastic-plastic transition strain
        epsX:            Threshold strain between compaction regimes
        kappa:           Compaction rate
        gammaS0:         Reference gamma at full density
        cS0:             Reference sound speed at full density
        c0:              Reference sound speed at initial porosity
        rhoS0:           Reference solid density
        jutziStateUpdate Apply state update rules from Jutzi 2008"""

    def pyinit1(self,
                nodeList = "const SolidNodeList<%(Dimension)s>&",
                phi0 = "const Field<%(Dimension)s, %(Dimension)s::Scalar>&",
                epsE = "const double",
                epsX = "const double",
                kappa = "const double",
                gammaS0 = "const double",
                cS0 = "const double",
                c0 = "const Field<%(Dimension)s, %(Dimension)s::Scalar>&",
                rhoS0 = "const double",
                jutziStateUpdate = "const bool"):
        """Constructor parameters:
        nodeList:        The NodeList we're going apply to
        phi0:            Initial porosity (field of values)
        epsE:            Elastic-plastic transition strain
        epsX:            Threshold strain between compaction regimes
        kappa:           Compaction rate
        gammaS0:         Reference gamma at full density
        cS0:             Reference sound speed at full density
        c0:              Reference sound speed at initial porosity
        rhoS0:           Reference solid density
        jutziStateUpdate Apply state update rules from Jutzi 2008"""

    #...........................................................................
    # Properties
    epsE = PYB11property(doc="Elastic-plastic transition strain")
    epsX = PYB11property(doc="Threshold strain between compaction regimes")
    kappa = PYB11property(doc="Compaction rate")
    gammaS0 = PYB11property(doc="Reference gamma at full density")
    strain = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    DstrainDt = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, StrainPorosity, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, StrainPorosity, virtual=True, pure_virtual=False)
