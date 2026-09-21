"""
Spheral ArtificialViscosity module.

Provides the artificial viscosity algorithms for use with the hydrodynamics methods.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"ArtificialViscosity/ArtificialViscosity.hh"',
                  '"ArtificialViscosity/ArtificialViscosityView.hh"',
                  '"ArtificialViscosity/MonaghanGingoldViscosity.hh"',
                  '"ArtificialViscosity/LimitedMonaghanGingoldViscosity.hh"',
                  '"ArtificialViscosity/MorrisMonaghanReducingViscosity.hh"',
                  '"ArtificialViscosity/CullenDehnenViscosity.hh"',
                  '"ArtificialViscosity/TensorMonaghanGingoldViscosity.hh"',
                  '"ArtificialViscosity/FiniteVolumeViscosity.hh"',
                  '"ArtificialViscosity/TensorSVPHViscosity.hh"',
                  '"ArtificialViscosity/TensorCRKSPHViscosity.hh"',
                  '"ArtificialViscosity/PythonArtificialViscosity.hh"',
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
from ArtificialViscosity import *
from MonaghanGingoldViscosity import *
from TensorMonaghanGingoldViscosity import *
from LimitedMonaghanGingoldViscosity import *
from MorrisMonaghanReducingViscosity import *
from CullenDehnenViscosity import *
from FiniteVolumeViscosity import *
from TensorSVPHViscosity import *
from TensorCRKSPHViscosity import *
from PythonArtificialViscosity import *

for ndim in dims:
    Dimension = f"Dim<{ndim}>"
    for pref in ["Artificial",
                 "MorrisMonaghanReducing",
                 "CullenDehnen",
                 "TensorSVPH",
                 "TensorCRKSPH",
                 "MonaghanGingold",
                 "TensorMonaghanGingold",
                 "LimitedMonaghanGingold",
                 "FiniteVolume"]:
        exec(f'''
{pref}Viscosity{ndim}d = PYB11TemplateClass({pref}Viscosity, template_parameters="{Dimension}")
''')

    exec(f'''
PythonScalarArtificialViscosity{ndim}d = PYB11TemplateClass(PythonArtificialViscosity, template_parameters=("{Dimension}", "{Dimension}::Scalar"))
PythonTensorArtificialViscosity{ndim}d = PYB11TemplateClass(PythonArtificialViscosity, template_parameters=("{Dimension}", "{Dimension}::Tensor"))
''')
