from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldListView
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralCompiledModules.SpheralFieldListView")
class FieldListView:

    PYB11typedefs = """
    using FieldListType = FieldList<%(Dimension)s, %(Value)s>;
    using FieldListViewType = FieldListView<%(Dimension)s, %(Value)s>;
    using FieldViewType = FieldView<%(Dimension)s, %(Value)s>;
    using Scalar = %(Dimension)s::Scalar;
"""

    def pyinitCopy(self, rhs="FieldListViewType&"):
        "Copy constructor"

    #...........................................................................
    # Methods
    @PYB11const
    def size(self):
        "Number of Fields"
        return "size_t"

    @PYB11const
    def empty(self):
        return "bool"

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
    @PYB11cppname("size")
    @PYB11const
    def __len__(self):
        return "size_t"

    @PYB11cppname("operator[]")
    @PYB11returnpolicy("reference")
    @PYB11keepalive(0,1)
    def __getitem__(self, index="size_t"):
        return "FieldViewType&"

    @PYB11returnpolicy("reference")
    @PYB11implementation("[](FieldListViewType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a FieldListView."

    @PYB11returnpolicy("reference")
    @PYB11const
    def __call__(self,
                 fieldIndex = "const size_t",
                 nodeIndex = "const size_t"):
        "Return the %(Value)s for the given (fieldIndex, nodeIndex)."
        return "%(Value)s&"

    #...........................................................................
    # Properties
    numFields = PYB11property("size_t", doc="Number of Fields")
    numElements = PYB11property("size_t", doc="Number of elements in all the associated Fields")
    numInternalElements = PYB11property("size_t", doc="Number of internal elements in all the associated Fields")
    numGhostElements = PYB11property("size_t", doc="Number of ghost elementes in all the associated Fields")
