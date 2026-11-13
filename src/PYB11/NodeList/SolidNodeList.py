from PYB11Generator import *
from FluidNodeList import FluidNodeList
from RestartMethods import *

#-------------------------------------------------------------------------------
# SolidNodeList template
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralCompiledModules.SpheralNodeList")
@PYB11dynamic_attr
class SolidNodeList(FluidNodeList):
    "Spheral SolidNodeList base class in %(Dimension)s, i.e.,  the NodeList for solid dynamics."

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef Field<%(Dimension)s, int> IntField;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
    typedef Field<%(Dimension)s, Vector> VectorField;
    typedef Field<%(Dimension)s, Tensor> TensorField;
    typedef Field<%(Dimension)s, SymTensor> SymTensorField;
"""

    def pyinit(self,
               name = "std::string",
               eos = "EquationOfState<%(Dimension)s>&",
               strength = "StrengthModel<%(Dimension)s>&",
               numInternal = ("size_t", "0u"),
               numGhost = ("size_t", "0u"),
               hmin = ("double", "1e-20"),
               hmax = ("double", "1e20"),
               hminratio = ("double", "0.1"),
               nPerh = ("double", "2.01"),
               maxNumNeighbors = ("size_t", "500u"),
               rhoMin = ("double", "1e-10"),
               rhoMax = ("double", "1e100")):
        "Constructor for a SolidNodeList class."
        return

    @PYB11virtual
    @PYB11const
    def soundSpeed(self, result="ScalarField&"):
        "Compute the current sound speed, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def bulkModulus(self, result="ScalarField&"):
        "Compute the current bulk modulus, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def shearModulus(self, result="ScalarField&"):
        "Compute the current shear modulus, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def yieldStrength(self, result="ScalarField&"):
        "Compute the current yield strength, storing the result in the argument ScalarField."
        return "void"

    @PYB11virtual
    @PYB11const
    def YoungsModulus(self,
                      result="ScalarField&",
                      K = "const ScalarField&",
                      mu = "const ScalarField&"):
        "Compute the current Youngs modulus, storing the result in the argument result, given the bulk modulus (K) and shear modulus (mu)."
        return "void"

    @PYB11virtual
    @PYB11const
    def longitudinalSoundSpeed(self,
                               result="ScalarField&",
                               rho = "const ScalarField&",
                               K = "const ScalarField&",
                               mu = "const ScalarField&"):
        """Compute the current longitudinal sound speed, storing the result in the argument result,
given the mass density, bulk modulus, and shear modulus (rho, K, mu)"""
        return "void"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def deviatoricStress(self):
        "The mass density field"
        return "const SymTensorField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def plasticStrain(self):
        "The plastic strain field"
        return "const ScalarField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def plasticStrainRate(self):
        "The plastic strain rate field"
        return "const ScalarField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def damage(self):
        "The damage field"
        return "const SymTensorField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def fragmentIDs(self):
        "The fragment IDs field"
        return "const IntField&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def particleTypes(self):
        "The particle type field"
        return "const IntField&"

    @PYB11implementation("[](SolidNodeList<%(Dimension)s>& self, const size_t nodeID) -> SymTensor& { return self.deviatoricStress()[nodeID]; }")
    @PYB11pycppname("deviatoricStress")
    @PYB11returnpolicy("reference_internal")
    def deviatoricStress1(self, nodeID="const size_t"):
        return "SymTensor&"

    @PYB11implementation("[](SolidNodeList<%(Dimension)s>& self, const size_t nodeID) -> Scalar& { return self.plasticStrain()[nodeID]; }")
    @PYB11pycppname("plasticStrain")
    @PYB11returnpolicy("reference_internal")
    def plasticStrain1(self, nodeID="const size_t"):
        return "Scalar&"

    @PYB11implementation("[](SolidNodeList<%(Dimension)s>& self, const size_t nodeID) -> Scalar& { return self.plasticStrainRate()[nodeID]; }")
    @PYB11pycppname("plasticStrainRate")
    @PYB11returnpolicy("reference_internal")
    def plasticStrainRate1(self, nodeID="const size_t"):
        return "Scalar&"

    @PYB11implementation("[](SolidNodeList<%(Dimension)s>& self, const size_t nodeID) -> SymTensor& { return self.damage()[nodeID]; }")
    @PYB11pycppname("damage")
    @PYB11returnpolicy("reference_internal")
    def damage1(self, nodeID="const size_t"):
        return "SymTensor&"

    @PYB11implementation("[](SolidNodeList<%(Dimension)s>& self, const size_t nodeID) -> int& { return self.fragmentIDs()[nodeID]; }")
    @PYB11returnpolicy("reference_internal")
    def fragmentID(self, nodeID="const size_t"):
        return "int&"

    @PYB11implementation("[](SolidNodeList<%(Dimension)s>& self, const size_t nodeID) -> int& { return self.particleTypes()[nodeID]; }")
    @PYB11returnpolicy("reference_internal")
    def particleType(self, nodeID="const size_t"):
        return "int&"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def strengthModel(self):
        "Return the strength model object this SolidNodeList is associated with."
        return "const StrengthModel<%(Dimension)s>&"

    # Comparison
    @PYB11pyname("__eq__")
    @PYB11implementation("[](const SolidNodeList<%(Dimension)s>& self, const SolidNodeList<%(Dimension)s>& rhs) { return static_cast<const NodeList<%(Dimension)s>&>(self) == static_cast<const NodeList<%(Dimension)s>&>(rhs); }")
    def equal(self, rhs="const SolidNodeList<%(Dimension)s>&"):
        "Equivalence test with another SolidNodeList"
        return "bool"

    @PYB11pyname("__ne__")
    @PYB11implementation("[](const SolidNodeList<%(Dimension)s>& self, const SolidNodeList<%(Dimension)s>& rhs) { return static_cast<const NodeList<%(Dimension)s>&>(self) != static_cast<const NodeList<%(Dimension)s>&>(rhs); }")
    def notEqual(self, rhs="const SolidNodeList<%(Dimension)s>&"):
        "Inequivalence test with another SolidNodeList"
        return "bool"

    @PYB11implementation("[](const SolidNodeList<%(Dimension)s>& self) -> std::uintptr_t { return reinterpret_cast<std::uintptr_t>(&self); }")
    def __hash__(self):
        "Make SolidNodeList objects hashable for Python"
        return "std::uintptr_t"

#-------------------------------------------------------------------------------
# Inject the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SolidNodeList)
