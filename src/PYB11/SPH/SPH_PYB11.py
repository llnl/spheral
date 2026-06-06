"""
Spheral SPH module.

Provides implementations of SPH, PSPH, and ASPH 
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from SPHBase import *
from SPH import *
from PSPH import *
from SolidSPH import *
from SolidSphericalSPH import *
from SPH_RAJA import *
from SolidSPH_RAJA import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"SPH/SPHBase.hh"',
                  '"SPH/SPH.hh"',
                  '"SPH/PSPH.hh"',
                  '"SPH/computeSPHSumMassDensity.hh"',
                  '"SPH/computeSPHOmegaGradhCorrection.hh"',
                  '"SPH/SPHRZ.hh"',
                  '"SPH/SphericalSPH.hh"',
                  '"SPH/SolidSPH.hh"',
                  '"SPH/SolidSPHRZ.hh"',
                  '"SPH/SolidSphericalSPH.hh"',
                  '"SPH/SPH_RAJA.hh"',
                  '"SPH/SolidSPH_RAJA.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"FileIO/FileIO.hh"',
                  '"ArtificialViscosity/ArtificialViscosity.hh"',
                  '"Neighbor/PairwiseField.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "KernelType")
def computeSPHSumMassDensity(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                             W = "const %(KernelType)s&",
                             sumOverAllNodeLists = "const bool",
                             position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                             mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                             H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                             massDensity = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the SPH mass density summation."
    return "void"

@PYB11template("Dimension")
def computeSPHOmegaGradhCorrection(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                                   W = "const TableKernel<%(Dimension)s>&",
                                   position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                                   H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                                   omegaGradh = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the SPH grad h correction due to Springel et al."
    return "void"

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    Dimension = f"Dim<{ndim}>"
    exec(f'''
SPHBase{ndim}d = PYB11TemplateClass(SPHBase, template_parameters="{Dimension}")
SPH{ndim}d = PYB11TemplateClass(SPH, template_parameters="{Dimension}")
PSPH{ndim}d = PYB11TemplateClass(PSPH, template_parameters="{Dimension}")
SolidSPH{ndim}d = PYB11TemplateClass(SolidSPH, template_parameters="{Dimension}")

SPH_RAJA{ndim}d = PYB11TemplateClass(SPH_RAJA, template_parameters="{Dimension}")
SolidSPH_RAJA{ndim}d = PYB11TemplateClass(SolidSPH_RAJA, template_parameters="{Dimension}")

computeSPHSumMassDensity{ndim}d = PYB11TemplateFunction(computeSPHSumMassDensity, template_parameters=("{Dimension}", "TableKernel<{Dimension}>"))
computeSPHOmegaGradhCorrection{ndim}d = PYB11TemplateFunction(computeSPHOmegaGradhCorrection, template_parameters="{Dimension}")
''')

if 1 in dims:
    from SphericalSPH import *
    computeSPHSumMassDensity1d_spherical = PYB11TemplateFunction(computeSPHSumMassDensity, template_parameters=("Dim<1>", "SphericalKernel"), pyname="computeSPHSumMassDensity1d")

if 2 in dims:
    from SPHRZ import *
    from SolidSPHRZ import *
