#-------------------------------------------------------------------------------
# State
#-------------------------------------------------------------------------------
from PYB11Generator import *
from StateBase import *

@PYB11template("Dimension")
class State(StateBase):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename StateBase<%(Dimension)s>::KeyType KeyType;
    typedef typename StateBase<%(Dimension)s>::FieldName FieldName;
    typedef typename StateBase<%(Dimension)s>::MeshPtr MeshPtr;
    typedef typename State<%(Dimension)s>::PackageList PackageList;
    typedef typename State<%(Dimension)s>::PolicyPointer PolicyPointer;
    using BoundaryType = Boundary<%(Dimension)s>;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                dataBase = "DataBase<%(Dimension)s>&",
                physicsPackages = "PackageList&"):
        "Construct using the Physics::registerState methods of the packages"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def enroll(self, field="FieldBase<%(Dimension)s>&"):
        "Enroll a field to track"
        return "void"

    # @PYB11virtual
    # @PYB11pycppname("enroll")
    # def enroll_sharedptr(self, field="std::shared_ptr<FieldBase<%(Dimension)s>>&"):
    #     "Enroll a shared_ptr<Field> to track"
    #     return "void"

    @PYB11virtual
    @PYB11pycppname("enroll")
    def enroll_fieldlist(self, fieldList="FieldListBase<%(Dimension)s>&"):
        "Enroll a FieldList to track"
        return "void"

    #...........................................................................
    # Operators
    def __eq__(self):
        return

    #...........................................................................
    # Methods
    def update(self,
               derivs = "StateDerivatives<%(Dimension)s>&",
               multiplier = "const double",
               t = "const double",
               dt = "const double",
               dependentOnly = "const bool"):
        "Advance using the derivatives object, assuming state1 = state0 + multiplier*derivs"
        return "void"

    @PYB11pycppname("enroll")
    def enroll1(self,
                key = "const KeyType&",
                policy = "PolicyPointer"):
        "Enroll an update policy under the given key"
        return "void"

    @PYB11pycppname("enroll")
    def enroll2(self,
                field = "FieldBase<%(Dimension)s>&",
                policy = "PolicyPointer"):
        "Enroll a Field and associated update policy"
        return "void"

    @PYB11pycppname("enroll")
    def enroll3(self,
                fieldList = "FieldListBase<%(Dimension)s>&",
                policy = "PolicyPointer"):
        "Enroll a FieldList and associated update policy"
        return "void"

    def removePolicy(self, key="const KeyType&"):
        "Remove the policy based on the lookup key"
        return "void"

    @PYB11pycppname("removePolicy")
    def removePolicy1(self, field="FieldBase<%(Dimension)s>&"):
        "Remove the policy associated with the Field"
        return "void"

    @PYB11pycppname("removePolicy")
    def removePolicy2(self,
                      fieldList="FieldListBase<%(Dimension)s>&",
                      clonePerField = "const bool"):
        "Remove the policy associated with the FieldList"
        return "void"

    @PYB11const
    @PYB11pycppname("policy")
    def policy1(self, key="const KeyType&"):
        "Get the update policy associated with a key"
        return "PolicyPointer"

    @PYB11const
    def serializeIndependentData(self,
                                 buf = "std::vector<double>&"):
        "Serialize any data associated with independent state"
        return "void"

    @PYB11const
    def deserializeIndependentData(self,
                                   buf = "const std::vector<double>&"):
        "Deserialize data associated with independent state"
        return "void"

    @PYB11const
    def serializeDerivatives(self,
                             buf = "std::vector<double>&",
                             derivs = "const StateDerivatives<%(Dimension)s>&"):
        "Serialize any derivatives associated with independent state"
        return "void"

    #...........................................................................
    # Template methods
    @PYB11template("Value")
    @PYB11const
    def policy(self, field="const Field<%(Dimension)s, %(Value)s>&"):
        "Return the update policy for the given Field"
        return "PolicyPointer"

    policy10 = PYB11TemplateMethod(policy, "int", pyname="policy")

    #...........................................................................
    # Proerties
    policyKeys = PYB11property("std::vector<KeyType>", "policyKeys", doc="The full set of keys for all policies")
    timeAdvanceOnly = PYB11property("bool", "timeAdvanceOnly", "timeAdvanceOnly", doc="Optionally trip a flag indicating policies should time advance only -- no replacing state!")
