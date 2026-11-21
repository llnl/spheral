"""
Spheral Boundary module.

Provides the Boundary base class and many boundary implementations.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from DistributedBoundary import *
from TreeDistributedBoundary import *
from BoundingVolumeDistributedBoundary import *
from RedistributeNodes import *
from SpaceFillingCurveRedistributeNodes import *
from MortonOrderRedistributeNodes import *
from PeanoHilbertOrderRedistributeNodes import *
from SortAndDivideRedistributeNodes import *
from VoronoiRedistributeNodes import *
from DistributeByXPosition import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Boundary/Boundary.hh"',
                  '"Distributed/DistributedBoundary.hh"',
                  '"Distributed/BoundingVolumeDistributedBoundary.hh"',
                  '"Distributed/TreeDistributedBoundary.hh"',
                  '"Distributed/RedistributeNodes.hh"',
                  '"Distributed/DistributeByXPosition.hh"',
                  '"Distributed/SpaceFillingCurveRedistributeNodes.hh"',
                  '"Distributed/MortonOrderRedistributeNodes.hh"',
                  '"Distributed/PeanoHilbertOrderRedistributeNodes.hh"',
                  '"Distributed/SortAndDivideRedistributeNodes1d.hh"',
                  '"Distributed/SortAndDivideRedistributeNodes2d.hh"',
                  '"Distributed/SortAndDivideRedistributeNodes3d.hh"',
                  '"Distributed/VoronoiRedistributeNodes.hh"',
                  '"Utilities/KeyTraits.hh"',
                  '"Field/Field.hh"',
                  '"Field/FieldList.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
DistributedBoundary%(ndim)id = PYB11TemplateClass(DistributedBoundary, template_parameters="%(Dimension)s")
TreeDistributedBoundary%(ndim)id = PYB11TemplateClass(TreeDistributedBoundary, template_parameters="%(Dimension)s")
BoundingVolumeDistributedBoundary%(ndim)id = PYB11TemplateClass(BoundingVolumeDistributedBoundary, template_parameters="%(Dimension)s")
RedistributeNodes%(ndim)id = PYB11TemplateClass(RedistributeNodes, template_parameters="%(Dimension)s")
SortAndDivideRedistributeNodesBase%(ndim)id = PYB11TemplateClass(SortAndDivideRedistributeNodes, template_parameters="%(Dimension)s")
SpaceFillngCurveRedistributeNodes%(ndim)id = PYB11TemplateClass(SpaceFillingCurveRedistributeNodes, template_parameters="%(Dimension)s")
MortonOrderRedistributeNodes%(ndim)id = PYB11TemplateClass(MortonOrderRedistributeNodes, template_parameters="%(Dimension)s")
PeanoHilbertOrderRedistributeNodes%(ndim)id = PYB11TemplateClass(PeanoHilbertOrderRedistributeNodes, template_parameters="%(Dimension)s")
VoronoiRedistributeNodes%(ndim)id = PYB11TemplateClass(VoronoiRedistributeNodes, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : ("Dim<" + str(ndim) +">")})

if 1 in dims:
    from SortAndDivideRedistributeNodes1d import *
    DistributeByXPosition1d = PYB11TemplateClass(DistributeByXPosition, template_parameters="Dim<1>")
    
if 2 in dims:
    from SortAndDivideRedistributeNodes2d import *
    DistributeByXPosition2d = PYB11TemplateClass(DistributeByXPosition, template_parameters="Dim<2>")
    
if 3 in dims:
    from SortAndDivideRedistributeNodes3d import *
