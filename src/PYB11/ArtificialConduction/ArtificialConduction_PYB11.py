"""
Spheral ArtificialConduction module.

Provides the artificial conduction methods.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from ArtificialConduction import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"ArtificialConduction/ArtificialConduction.hh"',
                  '"Field/FieldList.hh"',
                  '"DataBase/DataBase.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our dimensional types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
ArtificialConduction%(ndim)id = PYB11TemplateClass(ArtificialConduction, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
