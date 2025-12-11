from PYB11Generator import *
from PureReplaceBoundedState import *

@PYB11module("SpheralCompiledModules.SpheralDataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension", "ValueType")
class ReplaceBoundedState(PureReplaceBoundedState):

    PYB11typedefs = """
    using KeyType = typename ReplaceBoundedState<%(Dimension)s, %(ValueType)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self,
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
        "Build with no dependencies"
        return

    @PYB11implementation("[](py::list depends, const %(ValueType)s minValue, const %(ValueType)s maxValue) { auto result = make_policy<ReplaceBoundedState<%(Dimension)s, %(ValueType)s>>(minValue, maxValue); for (auto x: depends) result->addDependency(x.cast<std::string>()); return result; }")
    def pyinit1(self,
                depends = ("py::list", "py::list()"),
                minValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::lowest())"),
                maxValue = ("const %(ValueType)s", "%(ValueType)s(std::numeric_limits<double>::max())")):
        "Build with dependencies"
        return

    #...........................................................................
    # Methods
    @PYB11virtual
    def updateAsIncrement(self,
                          key = "const KeyType&",
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&",
                          multiplier = "const double",
                          t = "const double",
                          dt = "const double"):
        "Update a Field assoicated with the given key"
        return "void"
