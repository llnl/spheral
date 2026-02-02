#-------------------------------------------------------------------------------
# ThreePointDamagedNodeCoupling
#-------------------------------------------------------------------------------
from PYB11Generator import *
import NodeCoupling

@PYB11template("Dimension")
class ThreePointDamagedNodeCoupling(NodeCoupling.NodeCoupling):
    """A functor class encapsulating how we couple solid nodes in the presence of
multiple materials and damage.

This for uses the "three point" formalism, which allows damaged points to
cut communication between pairs that talk across them."""

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
"""

    def pyinit(self,
               state = "const State<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               pairs = "NodePairList&"):
        "Constructor"
