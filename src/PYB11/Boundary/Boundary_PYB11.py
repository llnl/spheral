"""
Spheral Boundary module.

Provides the Boundary base class and many boundary implementations.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Boundary/Boundary.hh"',
                  '"Boundary/PlanarBoundary.hh"',
                  '"Boundary/ReflectingBoundary.hh"',
                  '"Boundary/RigidBoundary.hh"',
                  '"Boundary/PeriodicBoundary.hh"',
                  '"Boundary/ConstantVelocityBoundary.hh"',
                  '"Boundary/ConstantXVelocityBoundary.hh"',
                  '"Boundary/ConstantYVelocityBoundary.hh"',
                  '"Boundary/ConstantZVelocityBoundary.hh"',
                  '"Boundary/ConstantRVelocityBoundary.hh"',
                  '"Boundary/ConstantBoundary.hh"',
                  '"Boundary/SphericalBoundary.hh"',
                  '"Boundary/CylindricalBoundary.hh"',
                  '"Boundary/AxisBoundaryRZ.hh"',
                  '"Boundary/SphericalOriginBoundary.hh"',
                  '"Boundary/CRKSPHVoidBoundary.hh"',
                  '"Boundary/InflowOutflowBoundary.hh"',
                  '"Boundary/mapPositionThroughPlanes.hh"',
                  '"Boundary/findNodesTouchingThroughPlanes.hh"',
                  '"Boundary/FacetedVolumeBoundary.hh"',
                  '"Field/Field.hh"',
                  '"Field/FieldList.hh"',
                  '"NodeList/NodeList.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
@PYB11template("Boundary1", "Boundary2")
@PYB11returnpolicy("reference")
@PYB11implementation("[](%(Boundary1)s* boundPtr) { return dynamic_cast<%(Boundary2)s*>(boundPtr); }")
def dynamicCastBoundary(boundPtr = "%(Boundary1)s*"):
    "Dynamic cast %(Boundary1)s* -> %(Boundary2)s*"
    return "%(Boundary2)s*"

@PYB11template("Dimension")
def mapPositionThroughPlanes(position = "const %(Dimension)s::Vector&",
                             enterPlane = "const GeomPlane<%(Dimension)s>&",
                             exitPlane = "const GeomPlane<%(Dimension)s>&"):
    """Function to take a position, and map it through the given enter/exit plane
pair and return the mapped position."""
    return "%(Dimension)s::Vector"

@PYB11template("Dimension")
def findNodesTouchingThroughPlanes(nodeList = "const NodeList<%(Dimension)s>&",
                                   enterPlane = "const GeomPlane<%(Dimension)s>&",
                                   exitPlane = "const GeomPlane<%(Dimension)s>&",
                                   hmultiplier = ("double", "1.0")):
    "Find the set of nodes that see through a pair of planes."
    return "std::vector<size_t>"

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from Boundary import *
from PlanarBoundary import *
from ReflectingBoundary import *
from RigidBoundary import *
from PeriodicBoundary import *
from ConstantVelocityBoundary import *
from ConstantXVelocityBoundary import *
from ConstantYVelocityBoundary import *
from ConstantZVelocityBoundary import *
from ConstantRVelocityBoundary import *
from ConstantBoundary import *
from CRKSPHVoidBoundary import *
from InflowOutflowBoundary import *
from FacetedVolumeBoundary import *

for ndim in dims:
    exec('''
Boundary%(ndim)id = PYB11TemplateClass(Boundary, template_parameters="%(Dimension)s")
PlanarBoundary%(ndim)id = PYB11TemplateClass(PlanarBoundary, template_parameters="%(Dimension)s")
ReflectingBoundary%(ndim)id = PYB11TemplateClass(ReflectingBoundary, template_parameters="%(Dimension)s")
RigidBoundary%(ndim)id = PYB11TemplateClass(RigidBoundary, template_parameters="%(Dimension)s")
PeriodicBoundary%(ndim)id = PYB11TemplateClass(PeriodicBoundary, template_parameters="%(Dimension)s")
ConstantVelocityBoundary%(ndim)id = PYB11TemplateClass(ConstantVelocityBoundary, template_parameters="%(Dimension)s")
ConstantXVelocityBoundary%(ndim)id = PYB11TemplateClass(ConstantXVelocityBoundary, template_parameters="%(Dimension)s")
ConstantBoundary%(ndim)id = PYB11TemplateClass(ConstantBoundary, template_parameters="%(Dimension)s")
CRKSPHVoidBoundary%(ndim)id = PYB11TemplateClass(CRKSPHVoidBoundary, template_parameters="%(Dimension)s")
InflowOutflowBoundary%(ndim)id = PYB11TemplateClass(InflowOutflowBoundary, template_parameters="%(Dimension)s")

vector_of_Boundary%(ndim)id = PYB11_bind_vector("Boundary<%(Dimension)s>*", opaque=True, local=False)

dynamicCastBoundary%(ndim)id = PYB11TemplateFunction(dynamicCastBoundary, template_parameters = ("Boundary<%(Dimension)s>", "PlanarBoundary<%(Dimension)s>"))

mapPositionThroughPlanes%(ndim)id = PYB11TemplateFunction(mapPositionThroughPlanes, template_parameters="%(Dimension)s", pyname="mapPositionThroughPlanes")
findNodesTouchingThroughPlanes%(ndim)id = PYB11TemplateFunction(findNodesTouchingThroughPlanes, template_parameters="%(Dimension)s", pyname="findNodesTouchingThroughPlanes")

''' % {"ndim"      : ndim,
       "Dimension" : ("Dim<" + str(ndim) +">")})

if 1 in dims:
    from SphericalOriginBoundary import *

if 2 in dims:
    ConstantYVelocityBoundary2d = PYB11TemplateClass(ConstantYVelocityBoundary, template_parameters="Dim<2>")
    FacetedVolumeBoundary2d = PYB11TemplateClass(FacetedVolumeBoundary, template_parameters="Dim<2>")
    from AxisBoundaryRZ import *

if 3 in dims:
    ConstantYVelocityBoundary3d = PYB11TemplateClass(ConstantYVelocityBoundary, template_parameters="Dim<3>")
    ConstantZVelocityBoundary3d = PYB11TemplateClass(ConstantZVelocityBoundary, template_parameters="Dim<3>")
    FacetedVolumeBoundary3d = PYB11TemplateClass(FacetedVolumeBoundary, template_parameters="Dim<3>")
    from CylindricalBoundary import *
    from SphericalBoundary import *

