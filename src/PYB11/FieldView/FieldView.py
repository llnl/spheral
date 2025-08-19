from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldView
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralFieldView")
class FieldView:

    PYB11typedefs = """
  using SelfType = FieldView<%(Dimension)s, %(Value)s>;
  using Scalar = typename %(Dimension)s::Scalar;
"""

    def pyinit(self, field="Field<%(Dimension)s, %(Value)s>&"):
        "Construct from a Field"

    #...........................................................................
    # Comparators
    def __eq__(self):
        return

    def __ne__(self):
        return

    @PYB11pycppname("__eq__")
    def __eq__S__(self, rhs="%(Value)s()"):
        "Equivalence comparision with a %(Value)s"
        return "bool"

    @PYB11pycppname("__ne__")
    def __ne__S__(self, rhs="%(Value)s()"):
        "Not equal comparision with a %(Value)s"
        return "bool"

    #...........................................................................
    # Sequence methods
    @PYB11cppname("numElements")
    @PYB11const
    def __len__(self):
        return "size_t"

    @PYB11cppname("operator[]")
    @PYB11returnpolicy("reference_internal")
    @PYB11implementation('[](SelfType& self, int i) { const auto n = self.numElements(); if (size_t(i) >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }')
    def __getitem__(self):
        #return "%(Value)s&"
        return

    @PYB11implementation("[](SelfType& self, int i, const %(Value)s v) { const auto n = self.numElements(); if (size_t(i) >= n) throw py::index_error(); self[(i %% n + n) %% n] = v; }")
    def __setitem__(self):
        "Set a value"

    @PYB11implementation("[](SelfType& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a FieldView."

    @PYB11returnpolicy("reference_internal")
    @PYB11implementation("[](SelfType& self, int i) { const auto n = self.numElements(); if (size_t(i) >= n) throw py::index_error(); return &self[(i %% n + n) %% n]; }")
    def __call__(self):
        "Index into a FieldView"
        return

    #...........................................................................
    # Properties
    numElements = PYB11property("size_t", doc="Number of elements in span")
    numInternalElements = PYB11property("size_t", doc="Number of internal elements in span")
    numGhostElements = PYB11property("size_t", doc="Number of ghost elements in span")
