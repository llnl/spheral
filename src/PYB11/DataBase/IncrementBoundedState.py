from PYB11Generator import *
from IncrementState import *

@PYB11module("SpheralCompiledModules.SpheralDataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension", "ValueType")
class IncrementBoundedState(IncrementState):

    PYB11typedefs = """
    using KeyType = typename IncrementBoundedState<%(Dimension)s, %(ValueType)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self,
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())"),
                wildCardDerivs = ("const bool", "false")):
        "Build with no dependencies"
        return

    @PYB11implementation("[](py::list depends, const %(ValueType)s minValue, const %(ValueType)s maxValue, const bool wildCardDerivs) { auto result = make_policy<IncrementBoundedState<%(Dimension)s, %(ValueType)s>>(minValue, maxValue, wildCardDerivs); for (auto x: depends) result->addDependency(x.cast<std::string>()); return result; }")
    def pyinit1(self,
                depends = ("py::list", "py::list()"),
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
        "Build with dependencies"
        return

    #...........................................................................
    # Methods
    @PYB11virtual
    def update(self,
               key = "const KeyType&",
               state = "State<%(Dimension)s>&",
               derivs = "StateDerivatives<%(Dimension)s>&",
               multiplier = "const double",
               t = "const double",
               dt = "const double"):
        "Update a Field assoicated with the given key"
        return "void"

    @PYB11static
    def prefix(self):
        "Prefix for key of derivatives"
        return "const std::string"

    @PYB11virtual
    @PYB11const
    def serializeData(self,
                      buf = "std::vector<double>&",
                      key = "const KeyType&",
                      state = "const State<%(Dimension)s>&"):
        "Serialize the data in the Field to a buffer"
        return "void"

    @PYB11virtual
    @PYB11const
    def deserializeData(self,
                        buf = "const std::vector<double>&",
                        key = "const KeyType&",
                        state = "const State<%(Dimension)s>&",
                        offset = "const size_t"):
        "Deserialize the data in the Field from a buffer"
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
    # Properties
    minValue = PYB11property(doc="Minimum bound for Field values")
    maxValue = PYB11property(doc="Maximum bound for Field values")
