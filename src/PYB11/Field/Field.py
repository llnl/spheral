from PYB11Generator import *
from FieldBase import FieldBase
from FieldView import FieldView

#-------------------------------------------------------------------------------
# Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralField")
class Field(FieldBase,
            FieldView):

    PYB11typedefs = """
    using SelfType = Field<%(Dimension)s, %(Value)s>;
    using ViewType = typename SelfType::ViewType;
"""

    def pyinit(self, name="std::string"):
        "Construct with a name"

    def pyinit1(self, name="std::string", field="const SelfType&"):
        "Construct as a copy of a Field with a new name"

    def pyinit2(self, name="std::string", nodeList="const NodeList<%(Dimension)s>&"):
        "Construct with a name and NodeList"

    def pyinit3(self,
                name = "std::string",
                nodeList = "const NodeList<%(Dimension)s>&",
                val = "%(Value)s"):
        "Construct with a name, NodeList, and initial value"

    # def pyinit4(self, rhs="const PYB11TrampolineField<%(Dimension)s, %(Value)s>&"):
    #     "Copy constructor"

    #...........................................................................
    # Comparators
    def __eq__(self):
        return

    def __ne__(self):
        return

    @PYB11pycppname("__eq__")
    def __eq___S__(self, rhs="%(Value)s()"):
        "Equivalence comparision with a %(Value)s"
        return "bool"

    @PYB11pycppname("__ne__")
    def __ne__S__(self, rhs="%(Value)s()"):
        "Not equal comparision with a %(Value)s"
        return "bool"

    # #...........................................................................
    # # Sequence methods
    # @PYB11cppname("size")
    # @PYB11const
    # def __len__(self):
    #     return "size_t"

    # @PYB11cppname("operator[]")
    # @PYB11returnpolicy("reference_internal")
    # @PYB11implementation('[](SelfType& self, int i) { const int n = self.size(); if (i >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }')
    # def __getitem__(self):
    #     #return "%(Value)s&"
    #     return

    # @PYB11implementation("[](SelfType& self, int i, const %(Value)s v) { const int n = self.size(); if (i >= n) throw py::index_error(); self[(i %% n + n) %% n] = v; }")
    # def __setitem__(self):
    #     "Set a value"

    # @PYB11implementation("[](const SelfType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
    # def __iter__(self):
    #     "Python iteration through a Field."

    # @PYB11returnpolicy("reference_internal")
    # @PYB11implementation("[](SelfType& self, int i) { const int n = self.size(); if (i >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }")
    # def __call__(self):
    #     "Index into a Field"
    #     #return "%(Value)s&"
    #     return

    #...........................................................................
    # FieldBase virtual methods
    @PYB11virtual
    @PYB11const
    def size(self):
        "Number of elements"
        return "size_t"

    @PYB11virtual
    def Zero(self):
        "Set all element values equal to zero"
        return "void"

    @PYB11virtual
    def setNodeList(self, nodeList="const NodeList<%(Dimension)s>&"):
        "Register this Field with the given NodeList"
        return "void"

    @PYB11virtual
    @PYB11const
    def packValues(self, nodeIDs="const std::vector<size_t>&"):
        "Serialize the indicated elements into a vector<char>"
        return "std::vector<char>"

    @PYB11virtual
    def unpackValues(self,
                     nodeIDs="const std::vector<size_t>&",
                     buffer = "const std::vector<char>&"):
        "Deserialize values from the given buffer"
        return "void"

    @PYB11virtual
    def copyElements(self,
                     fromIndices="const std::vector<size_t>&",
                     toIndices="const std::vector<size_t>&"):
        "Copy a range of values from/to elements of the Field"
        return "void"

    @PYB11virtual
    @PYB11const
    def fixedSizeDataType(self):
        return "bool"

    @PYB11virtual
    @PYB11const
    def numValsInDataType(self):
        return "size_t"

    @PYB11virtual
    @PYB11const
    def sizeofDataType(self):
        return "size_t"

    @PYB11virtual
    @PYB11const
    def computeCommBufferSize(packIndices = "const std::vector<size_t>&",
                              sendProc = "int",
                              recvProc = "int"):
        return "size_t"

    #...........................................................................
    # Methods
    @PYB11const
    @PYB11implementation("[](const Field<%(Dimension)s, %(Value)s>& self) -> py::bytes { auto buf = self.serialize(); return py::bytes(&(*buf.begin()), buf.size()); }")
    def serialize(self):
        "Serialize Field to a vector<char>/python bytes object"
        return "py::bytes"

    @PYB11implementation("""[](Field<%(Dimension)s, %(Value)s>& self, const py::bytes& byts) -> void {
        std::string bufstr(byts);
        std::vector<char> buf(bufstr.begin(), bufstr.end());
        self.deserialize(buf);
    }""")
    def deserialize(self, buf="py::bytes&"):
        "Deserialize from a vector<char>/python bytes object"
        return "void"

    @PYB11const
    @PYB11implementation("[](const Field<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.internalValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def internalValues(self):
        "Return a python list (as a copy) of just the internal values"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const Field<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.ghostValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def ghostValues(self):
        "Return a python list (as a copy) of just the ghost values"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const Field<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.allValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def allValues(self):
        "Return a python list (as a copy) of all values in the Field"
        return "py::list"

    def view(self):
        return "ViewType&"
