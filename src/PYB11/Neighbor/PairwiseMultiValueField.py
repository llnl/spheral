from PYB11Generator import *

#-------------------------------------------------------------------------------
# PairwiseField
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value", "size_t numElements")
@PYB11module("SpheralCompiledModules.SpheralNeighbor")
@PYB11cppname("PairwiseField")
class PairwiseMultiValueField:

    PYB11typedefs = """
    using SELF = PairwiseField<%(Dimension)s, %(Value)s, %(numElements)s>;
    using reference = typename SELF::reference;
"""

    def pyinit(self,
               connectivity = "const ConnectivityMap<%(Dimension)s>&"):
        "Constructor"

    def pyinit1(self,
                rhs = "const SELF&"):
        "Copy constructor"
        
    #...........................................................................
    # Sequence methods
    @PYB11cppname("size")
    @PYB11const
    def __len__(self):
        return "size_t"

    @PYB11cppname("operator[]")
    @PYB11returnpolicy("reference_internal")
    @PYB11implementation('[](SELF& self, int i) { const int n = self.size(); if (i >= n) throw py::index_error(); return self[(i %% n + n) %% n]; }')
    def __getitem__(self):
        return

    # @PYB11implementation("[](SELF& self, int i, const %(Value)s v) { const int n = self.size(); if (i >= n) throw py::index_error(); self[(i %% n + n) %% n] = v; }")
    # def __setitem__(self):
    #     "Set a value"

    @PYB11implementation("[](const SELF& self) { return py::make_iterator(self.begin(), self.end()); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a Field."

    @PYB11returnpolicy("reference_internal")
    @PYB11implementation("[](SELF& self, int i) { const int n = self.size(); if (i >= n) throw py::index_error(); return self[(i %% n + n) %% n]; }")
    def __call__(self):
        "Index into a Field"
        return

    #...........................................................................
    # Methods
    @PYB11returnpolicy("reference_internal")
    def __call__(self,
                 x = "const NodePairIdxType&"):
        "Index by NodePair"
        return "reference"

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def pairs(self):
        "NodePairList this PairwiseField is defined on"
        return "const NodePairList&"

    #...........................................................................
    # Properties
    size = PYB11property("size_t", doc="size of the PairwiseField")
