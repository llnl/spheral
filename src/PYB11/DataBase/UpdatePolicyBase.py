from PYB11Generator import *
from UpdatePolicyBase import *

@PYB11module("SpheralDataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension")
class UpdatePolicyBase:

    PYB11typedefs = """
    using KeyType = typename UpdatePolicyBase<%(Dimension)s>::KeyType;
"""

    # #...........................................................................
    # # Constructors
    # @PYB11implementation("[](py::list depends) { auto result = make_policy<UpdatePolicyBase<%(Dimension)s>>(); for (auto x: depends) result->addDependency(x.cast<std::string>()); return result; }")
    # def pyinit(self,
    #            depends = ("py::list", "py::list()")):
    #     "Build with the given dependencies"
    #     return

    #...........................................................................
    # Virtual methods
    @PYB11pure_virtual
    def update(self,
               key = "const KeyType&",
               state = "State<%(Dimension)s>&",
               derivs = "StateDerivatives<%(Dimension)s>&",
               multiplier = "const double",
               t = "const double",
               dt = "const double"):
        "Update a FieldList assoicated with the given key"
        return "void"

    @PYB11virtual
    def updateAsIncrement(self,
                          key = "const KeyType&",
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&",
                          multiplier = "const double",
                          t = "const double",
                          dt = "const double"):
        """An alternate method to be called when you want to specify that the derivative information
should be assumed to not necessarily be properly time-centered, and therefore you should 
only use time advancement ideas, no "replace" or more sophisticated approaches.
Default to just calling the generic method."""
        return "void"

    @PYB11virtual
    @PYB11const
    def clonePerField(self):
        "Returns whether this policy should be cloned for each Field in a FieldList or not"
        return "bool"

    @PYB11virtual
    @PYB11const
    def serializeData(self,
                      buf = "std::vector<double>&",
                      key = "const KeyType&",
                      state = "const State<%(Dimension)s>&"):
        "Serialize the data in the state to a buffer"
        return "void"

    @PYB11virtual
    @PYB11const
    def deserializeData(self,
                        buf = "const std::vector<double>&",
                        key = "const KeyType&",
                        state = "const State<%(Dimension)s>&",
                        offset = "const size_t"):
        "Deserialize the data in the state from a buffer"
        return "size_t"

    @PYB11virtual
    @PYB11const
    def serializeDerivatives(self,
                             buf = "std::vector<double>&",
                             key = "const KeyType&",
                             state = "const StateDerivatives<%(Dimension)s>&"):
        "Serialize the data in the derivatives to a buffer"
        return "void"

    #...........................................................................
    # Methods
    @PYB11virtual
    @PYB11const
    def independent(self):
        return "bool"

    @PYB11virtual
    @PYB11const
    def dependent(self):
        return "bool"

    @PYB11const
    def dependencies(self):
        "Return the set of field names that this state depends upon (if any)."
        return "const std::vector<std::string>&"

    def addDependency(self,
                      depend = "const std::string&"):
        "Allow the addition of new dependencies."
        return "void"

    @PYB11static
    def wildcard(self):
        "The wildcard string for comparing dependencies"
        return "std::string"

    #...........................................................................
    # Comparators
    @PYB11pure_virtual
    @PYB11const
    @PYB11cppname("operator==")
    @PYB11pyname("equals")      # Prevent converting this to a normal python comparison since it's a pure virtual function
    def __eq__(self,
               rhs = "const UpdatePolicyBase<%(Dimension)s>&"):
        return "bool"
