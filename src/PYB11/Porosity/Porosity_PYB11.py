"""
Spheral Porosity module.

Provides porosity models for modeling solid materials with unresolved void spaces
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from SolidEquationOfState import *

from PorosityModel import *
from StrainPorosity import *
from PalphaPorosity import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Porosity/PorosityModel.hh"',
                  '"Porosity/StrainPorosity.hh"',
                  '"Porosity/PalphaPorosity.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our dimensional types
#-------------------------------------------------------------------------------
for ndim in dims:
    Dimension = "Dim<" + str(ndim) + ">"
    exec(f'''
PorosityModel{ndim}d = PYB11TemplateClass(PorosityModel, template_parameters="{Dimension}")
StrainPorosity{ndim}d = PYB11TemplateClass(StrainPorosity, template_parameters="{Dimension}")
PalphaPorosity{ndim}d = PYB11TemplateClass(PalphaPorosity, template_parameters="{Dimension}")
''')
