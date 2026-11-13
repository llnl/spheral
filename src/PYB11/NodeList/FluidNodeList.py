from PYB11Generator import *
from NodeList import NodeList
from RestartMethods import *

#-------------------------------------------------------------------------------
# FluidNodeList template
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralNodeList")
@PYB11dynamic_attr
class FluidNodeList(NodeList):
    "Spheral FluidNodeList base class in %(Dimension)s, i.e.,  the NodeList for fluid hydrodynamics."

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
    typedef Field<%(Dimension)s, Vector> VectorField;
    typedef Field<%(Dimension)s, Tensor> TensorField;
    typedef Field<%(Dimension)s, SymTensor> SymTensorField;
"""

    def pyinit(self,
               name = "std::string",
               eos = "EquationOfState<%(Dimension)s>&",
               numInternal = ("size_t", "0u"),
               numGhost = ("size_t", "0u"),
               hmin = ("double", "1e-20"),
               hmax = ("double", "1e20"),
               hminratio = ("double", "0.1"),
               nPerh = ("double", "2.01"),
               maxNumNeighbors = ("size_t", "500u"),
               rhoMin = ("double", "1e-10"),
               rhoMax = ("double", "1e100")):
        "Constructor for a FluidNodeList class."
        return

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def massDensity(self):
        "The mass density field"
        return "const ScalarField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def massDensityRZ(self):
        "The mass density (RZ) field"
        return "const ScalarField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def specificThermalEnergy(self):
        "The specific thermal energy field"
        return "const ScalarField&"

    @PYB11implementation("[](FluidNodeList<%(Dimension)s>& self, const size_t nodeID) -> Scalar& { return self.massDensity()[nodeID]; }")
    @PYB11pycppname("massDensity")
    @PYB11returnpolicy("reference_internal")
    def massDensity1(self, nodeID="const size_t"):
        return "Scalar&"

    @PYB11implementation("[](FluidNodeList<%(Dimension)s>& self, const size_t nodeID) -> Scalar& { return self.massDensityRZ()[nodeID]; }")
    @PYB11pycppname("massDensityRZ")
    @PYB11returnpolicy("reference_internal")
    def massDensityRZ1(self, nodeID="const size_t"):
        return "Scalar&"

    @PYB11implementation("[](FluidNodeList<%(Dimension)s>& self, const size_t nodeID) -> Scalar& { return self.specificThermalEnergy()[nodeID]; }")
    @PYB11pycppname("specificThermalEnergy")
    @PYB11returnpolicy("reference_internal")
    def specificThermalEnergy1(self, nodeID="const size_t"):
        return "Scalar&"

    @PYB11pycppname("massDensity")
    def setmassDensity(self, val="const ScalarField&"):
        "Set the mass density field"
        return "void"

    @PYB11pycppname("massDensityRZ")
    def setmassDensityRZ(self, val="const ScalarField&"):
        "Set the mass density (RZ) field"
        return "void"

    @PYB11pycppname("specificThermalEnergy")
    def setspecificThermalEnergy(self, val="const ScalarField&"):
        "Set the specific thermal energy field"
        return "void"

    @PYB11virtual
    @PYB11const
    def pressure(self, result="ScalarField&"):
        "Compute the current pressure, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def temperature(self, result="ScalarField&"):
        "Compute the current temperature, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def soundSpeed(self, result="ScalarField&"):
        "Compute the current sound speed, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def volume(self, result="ScalarField&"):
        "Compute the current volume, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def linearMomentum(self, result="VectorField&"):
        "Compute the current linear momentum, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def totalEnergy(self, result="ScalarField&"):
        "Compute the current total energy, storing the result in the argument ScalarField."
        return "void"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def equationOfState(self):
        "Return the equation of state object this FluidNodeList is associated with."
        return "const EquationOfState<%(Dimension)s>&"

    @PYB11cppname("equationOfState")
    def setequationOfState(self, equationOfState="const EquationOfState<%(Dimension)s>&"):
        "Set the equation of state for this FluidNodeList."
        return "void"

    #...........................................................................
    # Comparison
    @PYB11pyname("__eq__")
    @PYB11implementation("[](const FluidNodeList<%(Dimension)s>& self, const FluidNodeList<%(Dimension)s>& rhs) { return static_cast<const NodeList<%(Dimension)s>&>(self) == static_cast<const NodeList<%(Dimension)s>&>(rhs); }")
    def equal(self, rhs="const FluidNodeList<%(Dimension)s>&"):
        "Equivalence test with another FluidNodeList"
        return "bool"

    @PYB11pyname("__ne__")
    @PYB11implementation("[](const FluidNodeList<%(Dimension)s>& self, const FluidNodeList<%(Dimension)s>& rhs) { return static_cast<const NodeList<%(Dimension)s>&>(self) != static_cast<const NodeList<%(Dimension)s>&>(rhs); }")
    def notEqual(self, rhs="const FluidNodeList<%(Dimension)s>&"):
        "Inequivalence test with another FluidNodeList"
        return "bool"

    @PYB11implementation("[](const FluidNodeList<%(Dimension)s>& self) -> std::uintptr_t { return reinterpret_cast<std::uintptr_t>(&self); }")
    def __hash__(self):
        "Make FluidNodeList objects hashable for Python"
        return "std::uintptr_t"

    #...........................................................................
    # Properties
    rhoMin = PYB11property("Scalar", "rhoMin", "rhoMin", doc="The minimum allowed mass density.")
    rhoMax = PYB11property("Scalar", "rhoMax", "rhoMax", doc="The maximum allowed mass density.")

#-------------------------------------------------------------------------------
# Inject the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, FluidNodeList)
