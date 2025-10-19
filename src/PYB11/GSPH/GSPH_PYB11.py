"""
Spheral GSPH module.

Provides implementations of Riemann SPH
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from GenericRiemannHydro import *
from GSPH import *
from MFM import *
from MFV import *
from WaveSpeeds import *
from Limiters import *
from RiemannSolvers import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"GSPH/GenericRiemannHydro.hh"',
                  '"GSPH/GSPH.hh"',
                  '"GSPH/MFM.hh"',
                  '"GSPH/MFV.hh"',
                  '"GSPH/WaveSpeeds/WaveSpeedBase.hh"',
                  '"GSPH/WaveSpeeds/AcousticWaveSpeed.hh"',
                  '"GSPH/WaveSpeeds/DavisWaveSpeed.hh"',
                  '"GSPH/WaveSpeeds/EinfeldtWaveSpeed.hh"',
                  '"GSPH/Limiters/LimiterBase.hh"',
                  '"GSPH/Limiters/MinModLimiter.hh"',
                  '"GSPH/Limiters/VanLeerLimiter.hh"',
                  '"GSPH/Limiters/VanAlbaLimiter.hh"',
                  '"GSPH/Limiters/SuperbeeLimiter.hh"',
                  '"GSPH/Limiters/OspreLimiter.hh"',
                  '"GSPH/Limiters/BarthJespersenLimiter.hh"',
                  '"GSPH/RiemannSolvers/RiemannSolverBase.hh"',
                  '"GSPH/RiemannSolvers/HLLC.hh"',
                  '"GSPH/RiemannSolvers/SecondOrderArtificialViscosity.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"Neighbor/PairwiseField.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
GradientType = PYB11enum(("RiemannGradient",
                          "HydroAccelerationGradient",
                          "SPHGradient",
                          "MixedMethodGradient",
                          "SPHSameTimeGradient",
                          "SPHUncorrectedGradient",
                          "NoGradient"), export_values = True)

NodeMotionType = PYB11enum(("Lagrangian",
                            "Eulerian",
                            "Fickian",
                            "XSPH",
                            "EulerianFickian"), export_values = False)

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
GenericRiemannHydro%(ndim)id = PYB11TemplateClass(GenericRiemannHydro, template_parameters="%(Dimension)s")
GSPH%(ndim)id = PYB11TemplateClass(GSPH, template_parameters="%(Dimension)s")
MFM%(ndim)id = PYB11TemplateClass(MFM, template_parameters="%(Dimension)s")
MFV%(ndim)id = PYB11TemplateClass(MFV, template_parameters="%(Dimension)s")
WaveSpeedBase%(ndim)id = PYB11TemplateClass(WaveSpeedBase, template_parameters="%(Dimension)s")
AcousticWaveSpeed%(ndim)id = PYB11TemplateClass(AcousticWaveSpeed, template_parameters="%(Dimension)s")
DavisWaveSpeed%(ndim)id = PYB11TemplateClass(DavisWaveSpeed, template_parameters="%(Dimension)s")
EinfeldtWaveSpeed%(ndim)id = PYB11TemplateClass(EinfeldtWaveSpeed, template_parameters="%(Dimension)s")
LimiterBase%(ndim)id = PYB11TemplateClass(LimiterBase, template_parameters="%(Dimension)s")
MinModLimiter%(ndim)id = PYB11TemplateClass(MinModLimiter, template_parameters="%(Dimension)s")
VanLeerLimiter%(ndim)id = PYB11TemplateClass(VanLeerLimiter, template_parameters="%(Dimension)s")
VanAlbaLimiter%(ndim)id = PYB11TemplateClass(VanAlbaLimiter, template_parameters="%(Dimension)s")
SuperbeeLimiter%(ndim)id = PYB11TemplateClass(SuperbeeLimiter, template_parameters="%(Dimension)s")
OspreLimiter%(ndim)id = PYB11TemplateClass(OspreLimiter, template_parameters="%(Dimension)s")
BarthJespersenLimiter%(ndim)id = PYB11TemplateClass(BarthJespersenLimiter, template_parameters="%(Dimension)s")
RiemannSolverBase%(ndim)id = PYB11TemplateClass(RiemannSolverBase, template_parameters="%(Dimension)s")
HLLC%(ndim)id = PYB11TemplateClass(HLLC, template_parameters="%(Dimension)s")
SecondOrderArtificialViscosity%(ndim)id = PYB11TemplateClass(SecondOrderArtificialViscosity, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

