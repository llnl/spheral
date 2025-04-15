"""
Spheral Gravity module.

Provides the implementations of gravitational physics packages in Spheral.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from Tree import *
from NBodyGravity import *
from TreeGravity import *
from PolyGravity import *
from ApproximatePolyhedralGravityModel import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Gravity/Tree.hh"',
                  '"Gravity/ApproximatePolyhedralGravityModel.hh"',
                  '"Physics/GenericBodyForce.hh"',
                  '"Gravity/NBodyGravity.hh"',
                  '"Gravity/TreeGravity.hh"',
                  '"Gravity/PolyGravity.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
GravityTimeStepType = PYB11enum(("AccelerationRatio", "DynamicalTime"), export_values=True,
                                doc="Algorithm choice for choosing timesteps for TreeGravity")

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
NBodyGravity%(ndim)id = PYB11TemplateClass(NBodyGravity, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

if 2 in dims:
    Tree2d = PYB11TemplateClass(Tree, template_parameters="Dim<2>")
    QuadTreeGravity = PYB11TemplateClass(TreeGravity, template_parameters="Dim<2>")
    PolyGravity2d = PYB11TemplateClass(PolyGravity, template_parameters="Dim<2>")

if 3 in dims:
    Tree3d = PYB11TemplateClass(Tree, template_parameters="Dim<3>")
    OctTreeGravity = PYB11TemplateClass(TreeGravity, template_parameters="Dim<3>")
    PolyGravity3d = PYB11TemplateClass(PolyGravity, template_parameters="Dim<3>")
