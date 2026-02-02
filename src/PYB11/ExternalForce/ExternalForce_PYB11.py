"""
Spheral ExternalForce module.

Provides the basis for force physics modules, and a few simple instances 
such as PointForce and the like.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Physics/GenericBodyForce.hh"',
                  '"ExternalForce/PointPotential.hh"',
                  '"ExternalForce/ConstantAcceleration.hh"',
                  '"ExternalForce/LinearAcceleration.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from PointPotential import *
from ConstantAcceleration import *
from LinearAcceleration import *

for ndim in dims:
    exec('''
PointPotential%(ndim)id = PYB11TemplateClass(PointPotential, template_parameters="%(Dimension)s")
ConstantAcceleration%(ndim)id = PYB11TemplateClass(ConstantAcceleration, template_parameters="%(Dimension)s")
LinearAcceleration%(ndim)id = PYB11TemplateClass(LinearAcceleration, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : ("Dim<" + str(ndim) +">")})
