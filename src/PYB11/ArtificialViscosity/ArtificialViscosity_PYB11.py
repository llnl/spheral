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
#from ArtificialViscosityView import *
from MonaghanGingoldViscosity import *
from TensorMonaghanGingoldViscosity import *
from LimitedMonaghanGingoldViscosity import *
from MorrisMonaghanReducingViscosity import *
from CullenDehnenViscosity import *
from FiniteVolumeViscosity import *
from TensorSVPHViscosity import *
from TensorCRKSPHViscosity import *

art_visc_names = ["MonaghanGingold", "TensorMonaghanGingold", "LimitedMonaghanGingold", "FiniteVolume"]

for ndim in dims:
    Dimension = f"Dim<{ndim}>"
    exec(f'''
ArtificialViscosity{ndim}d = PYB11TemplateClass(ArtificialViscosity, template_parameters="{Dimension}")
#ScalarArtificialViscosityView{ndim}d = PYB11TemplateClass(ArtificialViscosityView, template_parameters=("{Dimension}", "{Dimension}::Scalar"))
#TensorArtificialViscosityView{ndim}d = PYB11TemplateClass(ArtificialViscosityView, template_parameters=("{Dimension}", "{Dimension}::Tensor"))
MorrisMonaghanReducingViscosity{ndim}d = PYB11TemplateClass(MorrisMonaghanReducingViscosity, template_parameters="{Dimension}")
CullenDehnenViscosity{ndim}d = PYB11TemplateClass(CullenDehnenViscosity, template_parameters="{Dimension}")
TensorSVPHViscosity{ndim}d = PYB11TemplateClass(TensorSVPHViscosity, template_parameters="{Dimension}")
TensorCRKSPHViscosity{ndim}d = PYB11TemplateClass(TensorCRKSPHViscosity, template_parameters="{Dimension}")
''')

    for avn in art_visc_names:
        exec(f'''
{avn}Viscosity{ndim}d = PYB11TemplateClass({avn}Viscosity, template_parameters="{Dimension}")
#{avn}ViscosityView{ndim}d = PYB11TemplateClass({avn}ViscosityView, template_parameters="{Dimension}")
''')
