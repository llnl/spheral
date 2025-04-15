"""
Spheral DEM module.

Provides implementations of DEM
"""
from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from DEMBase import *
from LinearSpringDEM import *
from SolidBoundaries import *
from ContactIndex import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"DEM/DEMBase.hh"',
                  '"DEM/DEMFieldNames.hh"',
                  '"DEM/LinearSpringDEM.hh"',
                  '"DEM/ContactStorageLocation.hh"',
                  '"DEM/SolidBoundary/SolidBoundaryBase.hh"',
                  '"DEM/SolidBoundary/InfinitePlaneSolidBoundary.hh"',
                  '"DEM/SolidBoundary/RectangularPlaneSolidBoundary.hh"',
                  '"DEM/SolidBoundary/CircularPlaneSolidBoundary.hh"',
                  '"DEM/SolidBoundary/CylinderSolidBoundary.hh"',
                  '"DEM/SolidBoundary/SphereSolidBoundary.hh"',
                  '"DEM/SolidBoundary/ClippedSphereSolidBoundary.hh"',
                  '"DataOutput/registerWithRestart.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"Field/FieldList.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
DEMBase%(ndim)id = PYB11TemplateClass(DEMBase, template_parameters="%(Dimension)s")
LinearSpringDEM%(ndim)id = PYB11TemplateClass(LinearSpringDEM, template_parameters="%(Dimension)s")
SolidBoundaryBase%(ndim)id = PYB11TemplateClass(SolidBoundaryBase, template_parameters="%(Dimension)s")
InfinitePlaneSolidBoundary%(ndim)id = PYB11TemplateClass(InfinitePlaneSolidBoundary, template_parameters="%(Dimension)s")
RectangularPlaneSolidBoundary%(ndim)id = PYB11TemplateClass(RectangularPlaneSolidBoundary, template_parameters="%(Dimension)s")
CircularPlaneSolidBoundary%(ndim)id = PYB11TemplateClass(CircularPlaneSolidBoundary, template_parameters="%(Dimension)s")
CylinderSolidBoundary%(ndim)id = PYB11TemplateClass(CylinderSolidBoundary, template_parameters="%(Dimension)s")
SphereSolidBoundary%(ndim)id = PYB11TemplateClass(SphereSolidBoundary, template_parameters="%(Dimension)s")
ClippedSphereSolidBoundary%(ndim)id = PYB11TemplateClass(ClippedSphereSolidBoundary, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

#-------------------------------------------------------------------------------
# expose our field names
#-------------------------------------------------------------------------------
class DEMFieldNames:
    particleRadius = PYB11readonly(static=True, returnpolicy="copy")
    momentOfInertia = PYB11readonly(static=True, returnpolicy="copy")
    compositeParticleIndex = PYB11readonly(static=True, returnpolicy="copy")
    angularVelocity = PYB11readonly(static=True, returnpolicy="copy")
    uniqueIndices = PYB11readonly(static=True, returnpolicy="copy")
    isActiveContact = PYB11readonly(static=True, returnpolicy="copy")
    neighborIndices = PYB11readonly(static=True, returnpolicy="copy")
    shearDisplacement = PYB11readonly(static=True, returnpolicy="copy")
    rollingDisplacement = PYB11readonly(static=True, returnpolicy="copy")
    torsionalDisplacement = PYB11readonly(static=True, returnpolicy="copy")
    equilibriumOverlap = PYB11readonly(static=True, returnpolicy="copy")
    solidBoundaries = PYB11readonly(static=True, returnpolicy="copy")
