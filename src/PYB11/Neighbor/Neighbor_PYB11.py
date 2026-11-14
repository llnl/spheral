"""
Spheral Neighbor module.

Provides the base class and implementations for neighbor finding in Spheral.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()
import SpheralConfigs

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/GeomPlane.hh"',
                  '"Boundary/Boundary.hh"',
                  '"Neighbor/Neighbor.hh"',
                  '"Neighbor/TreeNeighbor.hh"',
                  '"Neighbor/ConnectivityMap.hh"',
                  '"Neighbor/NodePairIdxType.hh"',
                  '"Neighbor/NodePairList.hh"',
                  '"Neighbor/PairwiseField.hh"',
                  '<utility>']

if not SpheralConfigs.gpu_enabled():
    PYB11includes += ['"Neighbor/PythonNeighbor.hh"']
                  

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
from Neighbor import *
from TreeNeighbor import *
from ConnectivityMap import *
from NodePairIdxType import *
from NodePairList import *
from PairwiseField import *
from PairwiseMultiValueField import *

if not SpheralConfigs.gpu_enabled():
    from PythonNeighbor import *

for ndim in dims:
    suffix = f"{ndim}d"
    Dimension = f"Dim<{ndim}>"
    exec(f'''
Neighbor{suffix} = PYB11TemplateClass(Neighbor, template_parameters="{Dimension}")
TreeNeighbor{suffix} = PYB11TemplateClass(TreeNeighbor, template_parameters="{Dimension}")
ConnectivityMap{suffix} = PYB11TemplateClass(ConnectivityMap, template_parameters="{Dimension}")
''')

    if not SpheralConfigs.gpu_enabled():
        exec(f'''
PythonNeighbor{suffix} = PYB11TemplateClass(PythonNeighbor, template_parameters="{Dimension}")
''')

    for Value in ("Scalar", "Vector", "Tensor", "SymTensor"):
        V = f"{Dimension}::{Value}"
        exec(f'''
{Value}PairwiseField{suffix} = PYB11TemplateClass(PairwiseField, template_parameters=("{Dimension}", "{V}", "1u"))
''')

    exec(f'''
PairwiseFieldScalarScalar{suffix} = PYB11TemplateClass(PairwiseMultiValueField, template_parameters=("{Dimension}", "{Dimension}::Scalar", "2u"))
''')
