"""
Spheral FSISPH module.

"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()


from SolidFSISPH import *
from SlideSurface import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"FSISPH/SolidFSISPH.hh"',
                  '"FSISPH/FSIFieldNames.hh"',
                  '"FSISPH/SlideSurface.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"FileIO/FileIO.hh"',
                  '"ArtificialViscosity/ArtificialViscosity.hh"',
                  '"Neighbor/PairwiseField.hh"']

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
InterfaceMethod = PYB11enum(("HLLCInterface", 
                             "ModulusInterface",
                             "NoInterface"), export_values = True)

KernelAveragingMethod = PYB11enum(("NeverAverageKernels", 
                                   "AlwaysAverageKernels",
                                   "AverageInterfaceKernels"), export_values = True)

FSIMassDensityMethod = PYB11enum(("FSISumMassDensity",
                                  "PressureCorrectSumMassDensity",
                                  "HWeightedSumMassDensity"), export_values = True)

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
SlideSurface%(ndim)id = PYB11TemplateClass(SlideSurface, template_parameters="%(Dimension)s")
SolidFSISPH%(ndim)id = PYB11TemplateClass(SolidFSISPH, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})


#-------------------------------------------------------------------------------
# expose our field names
#-------------------------------------------------------------------------------
class FSIFieldNames:
    pressureGradient = PYB11readonly(static=True, returnpolicy="copy")
    specificThermalEnergyGradient = PYB11readonly(static=True, returnpolicy="copy")
    interfaceFlags = PYB11readonly(static=True, returnpolicy="copy")
    interfaceAreaVectors = PYB11readonly(static=True, returnpolicy="copy")
    interfaceNormals = PYB11readonly(static=True, returnpolicy="copy")
    interfaceAngles = PYB11readonly(static=True, returnpolicy="copy")
    interfaceFraction = PYB11readonly(static=True, returnpolicy="copy")
    interfaceSmoothness = PYB11readonly(static=True, returnpolicy="copy")
    smoothedInterfaceNormals = PYB11readonly(static=True, returnpolicy="copy")
    interfaceSmoothnessNormalization = PYB11readonly(static=True, returnpolicy="copy")
