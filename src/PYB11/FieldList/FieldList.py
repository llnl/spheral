from PYB11Generator import *
from FieldListBase import *
from FieldListView import FieldListView

#-------------------------------------------------------------------------------
# FieldList
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralFieldList")
class FieldList(FieldListBase,
                FieldListView):

    PYB11typedefs = """
    using FieldListType = FieldList<%(Dimension)s, %(Value)s>;
    using FieldType = Field<%(Dimension)s, %(Value)s>;
    using FieldListViewType = FieldListView<%(Dimension)s, %(Value)s>;
    using FieldViewType = FieldView<%(Dimension)s, %(Value)s>;
    using NodeListType = NodeList<%(Dimension)s>;
    using Vector = %(Dimension)s::Vector;
    using SymTensor = %(Dimension)s::SymTensor;
    using ViewType = typename FieldListType::ViewType;
"""

    def pyinit(self):
        "Default constructor"

    def pyinit1(self, name="FieldStorageType"):
        "Build with the storage type for the fields"

    def pyinitCopy(self, rhs="const FieldListType&"):
        "Copy constructor"

    #...........................................................................
    # Methods
    def storageType(self):
        "The method whereby Fields are stored/referenced by this FieldList"
        return

    def copyFields(self):
        """Force the the FieldList to make copies of all the Fields it currently\\npoints to and point to the internally cached copies instead.\\nAlso sets storageType to CopyFields."""
        return "void"

    @PYB11pycppname("copyFields")
    def copyFieldsF(self, fieldList="const FieldListType&"):
        "Store copies of Fields from another FieldList"
        return "void"

    @PYB11const
    def haveField(self, field="const FieldType&"):
        "Test if the specified field is pointed to by this FieldList"
        return "bool"

    @PYB11const
    def haveNodeList(self, nodeList="const NodeListType&"):
        "Test if there are any Fields based on the given NodeList currently pointed\\nto by this FieldList"
        return "bool"

    def assignFields(self, fieldList="const FieldListType&"):
        "Set all Fields pointed to by this FieldList equal to those of the given FieldList."
        return "void"

    def referenceFields(self, fieldList="const FieldListType&"):
        "Make this FieldList reference the Fields of another (no copying of Field data)."
        return "void"

    def appendField(self, field="const FieldType&"):
        "Add the given Field to this FieldList.\\nstorageType==ReferenceFields -> just make a new member Field pointing to the given Field.\\nstorageType==CopyFields -> copy the Field into a locally cached\\nversion, and add a\\n                           reference pointing to that local copy"
        return "void"

    def deleteField(self, field="const FieldType&"):
        "Remove the given Field from this FieldList."
        return "void"

    def appendNewField(self,
                       name = "const std::string",
                       nodeList = "const NodeListType&",
                       value = "const %(Value)s"):
        "Create a new field with the given (name, value) and add it to the FieldList.\\nNote this only makes sense when storageType==CopyFields"
        return "void"

    @PYB11returnpolicy("reference_internal")
    @PYB11implementation("[](const FieldListType& self, const NodeListType& nodeList) { return *self.fieldForNodeList(nodeList); }")
    def fieldForNodeList(self, nodeList="const NodeListType&"):
        "Return the Field associated with the given NodeList"
        return "FieldType*"

    @PYB11const
    def setMasterNodeLists(self,
                           r = "const Vector&",
                           H = "const SymTensor&",
                           masterLists = "std::vector<std::vector<int>>&",
                           coarseNeighbors = "std::vector<std::vector<int>>&"):
        "Set the master neighbor information based on (r,H)"
        return "void"

    @PYB11const
    def setMasterNodeLists(self,
                           r = "const Vector&",
                           masterLists = "std::vector<std::vector<int>>&",
                           coarseNeighbors = "std::vector<std::vector<int>>&"):
        "Set the master neighbor information based on position r (assume zero associated length scale)"
        return "void"

    @PYB11const
    def setRefineNodeLists(self,
                           r = "const Vector&",
                           H = "const SymTensor&",
                           coarseNeighbors = "const std::vector<std::vector<int>>&",
                           refineNeighbors = "std::vector<std::vector<int>>&"):
        "Set the refine neighbor information based on (r,H)"
        return "void"

    @PYB11const
    def setRefineNodeLists(self,
                           r = "const Vector&",
                           coarseNeighbors = "const std::vector<std::vector<int>>&",
                           refineNeighbors = "std::vector<std::vector<int>>&"):
        "Set the refine neighbor information based on r (assume zero associated length scale)"
        return "void"

    def Zero(self):
        "Set all element values equal to zero"
        return "void"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def nodeListPtrs(self):
        "The NodeLists for Fields in this FieldList"
        return "const std::vector<NodeListType*>&"

    @PYB11const
    @PYB11implementation("[](const FieldList<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.internalValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def internalValues(self):
        "Return a python list (as a copy) of just the internal values"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const FieldList<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.ghostValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def ghostValues(self):
        "Return a python list (as a copy) of just the ghost values"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const FieldList<%(Dimension)s, %(Value)s>& self) -> py::list { const auto vals = self.allValues(); py::list result; for (const auto& x: vals) result.append(x); return result; }")
    def allValues(self):
        "Return a python list (as a copy) of all values in the FieldList"
        return "py::list"

    def view(self):
        return "ViewType"

    #...........................................................................
    # Comparators
    def __eq__(self):
        return

    def __ne__(self):
        return

    def __eq__(self, rhs="%(Value)s()"):
        "Equivalence comparision with a %(Value)s"
        return "bool"

    def __ne__(self, rhs="%(Value)s()"):
        "Not equal comparision with a %(Value)s"
        return "bool"

    #...........................................................................
    # Sequence methods
    # @PYB11cppname("size")
    # @PYB11const
    # def __len__(self):
    #     return "size_t"

    @PYB11cppname("operator[]")
    @PYB11returnpolicy("reference")
    @PYB11keepalive(0,1)
    def __getitem__(self, index="const size_t"):
        return "FieldType*"

    @PYB11returnpolicy("reference")
    @PYB11implementation("[](const FieldListType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a FieldList."

    # def __call__(self,
    #              fieldIndex = "const size_t",
    #              nodeIndex = "const size_t"):
    #     "Return the %(Value)s for the given (fieldIndex, nodeIndex)."
    #     return "%(Value)s&"

    #...........................................................................
    # Properties
    storageType = PYB11property("FieldStorageType", doc="The method whereby Fields are stored/referenced by this FieldList")
