#-------------------------------------------------------------------------------
# PalphaPorosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from PorosityModel import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralPorosity")
class PalphaPorosity(PorosityModel):
    """An implementation of the P-alpha porosity model described in

Jutzi, M., Benz, W., & Michel, P. (2008). Numerical simulations of impacts
involving porous bodies.I. Implementing sub-resolution porosity in a 3D SPH
hydrocode. Icarus, 198(1), 242–255.

This model assumes you will provide a solid EOS which will be modified.
The underlying actualy solid EOS should provide the reference density, which
will be treated here as the compacted true solid reference density.

Note this model introduces a new state variable, the distention (alpha), which
the pressure now depends on.  This implies our usual definition of P(rho, eps)
now becomes P(rho, eps, alpha).  Our EOS interface does not recognize this
this parameter, so we store alpha locally and only allow Field updates of the
pressure (forbidding the single value P lookup the EOS usually allows).

PalphaPorosity is the physics module which time evolves the distention 
parameter (alpha) and gives it to the PorousEquationOfState.
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
               Pe = "const double",
               Pt = "const double",
               Ps = "const double",
               alphae = "const double",
               alphat = "const double",
               n1 = "const double",
               n2 = "const double",
               cS0 = "const double",
               c0 = "const double",
               rhoS0 = "const double",
               jutziStateUpdate = "const bool"):
        """Constructor parameters:
        nodeList:         The SolidNodeList we're going apply to
        phi0:             Initial porosity (single value)
        Pe:               Elastic pressure threshold
        Pt:               Transition pressure (Pe <= Pt)
        Ps:               Solid transition pressure (from fit, Pt <= Ps)
        alphae:           alpha at P=Pe
        alphat:           Transition distension
        n1:               Fitted exponent for plastic distention evolution
        n2:               Fitted exponent for plastic distention evolution
        cS0:              Reference sound speed at full density
        c0:               Reference sound speed at initial porosity
        rhoS0:            Reference solid density
        jutziStateUpdate: Optionally update the deviatoric stress and damage as described in Jutzi 2008"""

    def pyinit1(self,
                nodeList = "const SolidNodeList<%(Dimension)s>&",
                phi0 = "const Field<%(Dimension)s, %(Dimension)s::Scalar>&",
                Pe = "const double",
                Pt = "const double",
                Ps = "const double",
                alphae = "const double",
                alphat = "const double",
                n1 = "const double",
                n2 = "const double",
                cS0 = "const double",
                c0 = "const Field<%(Dimension)s, %(Dimension)s::Scalar>&",
                rhoS0 = "const double",
                jutziStateUpdate = "const bool"):
        """Constructor parameters:
        nodeList:         The SolidNodeList we're going apply to
        phi0:             Initial porosity (field of values)
        Pe:               Elastic pressure threshold
        Pt:               Transition pressure (Pe <= Pt)
        Ps:               Solid transition pressure (from fit, Pt <= Ps)
        alphae:           alpha at P=Pe
        alphat:           Transition distension
        n1:               Fitted exponent for plastic distention evolution
        n2:               Fitted exponent for plastic distention evolution
        cS0:              Reference sound speed at full density
        c0:               Reference sound speed at initial porosity (field)
        rhoS0:            Reference solid density
        jutziStateUpdate: Optionally update the deviatoric stress and damage as described in Jutzi 2008"""

    #...........................................................................
    # Methods
    @PYB11const
    def phi(self):
        "Compute the current porosity from the distention"
        return "ScalarField"

    #...........................................................................
    # Properties
    Pe = PYB11property(doc="Elastic pressure threshold")
    Pt = PYB11property(doc="Transition pressure (Pe <= Pt)")
    Ps = PYB11property(doc="Transition pressure (Pe <= Pt)")
    alphae = PYB11property(doc="Elastic distension threshold")
    alphat = PYB11property(doc="Elastic distension threshold")
    n1 = PYB11property(doc="Fitted exponent for plastic distention evolution")
    n2 = PYB11property(doc="Fitted exponent for plastic distention evolution")
    partialPpartialEps = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    partialPpartialRho = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, PalphaPorosity, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, PalphaPorosity, virtual=True, pure_virtual=False)
