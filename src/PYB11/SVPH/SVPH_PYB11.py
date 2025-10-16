"""
Spheral SVPH module.

Provides implementations of SVPH (Smoothed Voronoi Particle Hydrodynamics)
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from SVPHFieldNames import *
from SVPHFacetedHydroBase import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"SVPH/SVPHFieldNames.hh"',
                  '"SVPH/sampleFieldListSVPH.hh"',
                  '"SVPH/gradientFieldListSVPH.hh"',
                  '"SVPH/SVPHFacetedHydroBase.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"Neighbor/ConnectivityMap.hh"',
                  '"FileIO/FileIO.hh"',
                  '"ArtificialViscosity/ArtificialViscosity.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "DataType")
def sampleFieldListSVPH(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
                        position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                        Hfield = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                        connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                        W = "const TableKernel<%(Dimension)s>&",
                        mesh = "const Mesh<%(Dimension)s>&",
                        firstOrderConsistent = "const bool"):
    "Use SVPH to sample a FieldList."
    return "FieldList<%(Dimension)s, %(DataType)s>"

@PYB11template("Dimension", "DataType")
def gradientFieldListSVPH(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
                          position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                          Hfield = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                          connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                          W = "const TableKernel<%(Dimension)s>&",
                          mesh = "const Mesh<%(Dimension)s>&",
                          firstOrderConsistent = "const bool"):
    "Use SVPH to take the gradient of a FieldList."
    return "FieldList<%(Dimension)s, typename MathTraits<%(Dimension)s, %(DataType)s>::GradientType>"

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    Dimension = f"Dim<{ndim}>"
    exec(f'''
SVPHFacetedHydroBase{ndim}d = PYB11TemplateClass(SVPHFacetedHydroBase, template_parameters="{Dimension}")
''')

    # SVPH interpolation
    for element in (f"Dim<{ndim}>::Scalar",
                    f"Dim<{ndim}>::Vector",
                    f"Dim<{ndim}>::Tensor",
                    f"Dim<{ndim}>::SymTensor"):
        label = PYB11mangle(element)
        exec(f'''
sampleFieldListSVPH{label} = PYB11TemplateFunction(sampleFieldListSVPH, template_parameters=("{Dimension}", "{element}"), pyname="sampleFieldListSVPH")
''')

    # SVPH gradient
    for element in (f"Dim<{ndim}>::Scalar",
                    f"Dim<{ndim}>::Vector"):
        label = PYB11mangle(element)
        exec(f'''
gradientFieldListSVPH{label} = PYB11TemplateFunction(gradientFieldListSVPH, template_parameters=("{Dimension}", "{element}"), pyname="gradientFieldListSVPH")
''')

