"""
Spheral Neighbor module.

Provides the base class and implementations for neighbor finding in Spheral.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/GeomPlane.hh"',
                  '"Boundary/Boundary.hh"',
                  '"Neighbor/GridCellIndex.hh"',
                  '"Neighbor/GridCellPlane.hh"',
                  '"Neighbor/Neighbor.hh"',
                  '"Neighbor/NestedGridNeighbor.hh"',
                  '"Neighbor/TreeNeighbor.hh"',
                  '"Neighbor/ConnectivityMap.hh"',
                  '"Neighbor/NodePairIdxType.hh"',
                  '"Neighbor/NodePairList.hh"',
                  '"Neighbor/PairwiseField.hh"',
                  '<utility>']
                  

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
NeighborSearchType = PYB11enum(("Gather", "Scatter", "GatherScatter"), export_values=True,
                               doc="Choice for how nodes should look for neighbors")

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from GridCellIndex import *
from GridCellPlane import *
from Neighbor import *
from NestedGridNeighbor import *
from TreeNeighbor import *
from ConnectivityMap import *
from NodePairIdxType import *
from NodePairList import *
from PairwiseField import *
from PairwiseMultiValueField import *

for ndim in dims:
    suffix = f"{ndim}d"
    Dimension = f"Dim<{ndim}>"
    exec(f'''
GridCellIndex{suffix} = PYB11TemplateClass(GridCellIndex, template_parameters="{Dimension}")
GridCellPlane{suffix} = PYB11TemplateClass(GridCellPlane, template_parameters="{Dimension}")
Neighbor{suffix} = PYB11TemplateClass(Neighbor, template_parameters="{Dimension}")
NestedGridNeighbor{suffix} = PYB11TemplateClass(NestedGridNeighbor, template_parameters="{Dimension}")
TreeNeighbor{suffix} = PYB11TemplateClass(TreeNeighbor, template_parameters="{Dimension}")
ConnectivityMap{suffix} = PYB11TemplateClass(ConnectivityMap, template_parameters="{Dimension}")

vector_of_GridCellIndex{suffix} = PYB11_bind_vector("GridCellIndex<{Dimension}>", opaque=True, local=False)
vector_of_vector_of_GridCellIndex{suffix} = PYB11_bind_vector("std::vector<GridCellIndex<{Dimension}>>", opaque=True, local=False)
''')

    for Value in ("Scalar", "Vector", "Tensor", "SymTensor"):
        V = f"{Dimension}::{Value}"
        exec(f'''
{Value}PairwiseField{suffix} = PYB11TemplateClass(PairwiseField, template_parameters=("{Dimension}", "{V}", "1u"))
''')

    exec(f'''
PairwiseFieldScalarScalar{suffix} = PYB11TemplateClass(PairwiseMultiValueField, template_parameters=("{Dimension}", "{Dimension}::Scalar", "2u"))
''')
