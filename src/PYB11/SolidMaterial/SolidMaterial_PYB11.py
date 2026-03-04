"""
Spheral SolidMaterial module.

Provides equations of state and material models for solids in Spheral.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from SolidEquationOfState import *

from LinearPolynomialEquationOfState import *
from GruneisenEquationOfState import *
from OsborneEquationOfState import *
from TillotsonEquationOfState import *
from MurnaghanEquationOfState import *
from SingularityEquationOfState import *
from SingularityEquationOfStateModels import *

from StrengthModel import *
from NullStrength import *
from ConstantStrength import *
from SteinbergGuinanStrength import *
from JohnsonCookStrength import *
from CollinsStrength import *
from iSALEROCKStrength import *

from PhysicsEvolvingMaterialLibrary import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"SolidMaterial/SolidEquationOfState.hh"',
                  '"SolidMaterial/LinearPolynomialEquationOfState.hh"',
                  '"SolidMaterial/GruneisenEquationOfState.hh"',
                  '"SolidMaterial/OsborneEquationOfState.hh"',
                  '"SolidMaterial/TillotsonEquationOfState.hh"',
                  '"SolidMaterial/MurnaghanEquationOfState.hh"',
                  '"SolidMaterial/StrengthModel.hh"',
                  '"SolidMaterial/ConstantStrength.hh"',
                  '"SolidMaterial/NullStrength.hh"',
                  '"SolidMaterial/PolynomialFit.hh"',
                  '"SolidMaterial/SteinbergGuinanStrength.hh"',
                  '"SolidMaterial/SteinbergGuinanLundStrength.hh"',
                  '"SolidMaterial/JohnsonCookStrength.hh"',
                  '"SolidMaterial/CollinsStrength.hh"',
                  '"SolidMaterial/iSALEROCKStrength.hh"',
                  '"SolidMaterial/PhysicsEvolvingMaterialLibrary.hh"',
                  '"SolidMaterial/SingularityEquationOfState.hh"',
                  '"SolidMaterial/SingularityEquationOfStateModels.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# NinthOrderPolynomialFit
#-------------------------------------------------------------------------------
class NinthOrderPolynomialFit:
    "Used by Steinberg-Guinnan strength to approximate cold and melt energies"

    def pyinit(self,
               C0 = "const double",
               C1 = "const double",
               C2 = "const double",
               C3 = "const double",
               C4 = "const double",
               C5 = "const double",
               C6 = "const double",
               C7 = "const double",
               C8 = "const double",
               C9 = "const double"):
        "Construct with coefficients"

    @PYB11const
    def __call__(self, x="const double"):
        return "double"

#-------------------------------------------------------------------------------
# Instantiate our dimensional types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
SolidEquationOfState%(ndim)id = PYB11TemplateClass(SolidEquationOfState, template_parameters="%(Dimension)s")
StrengthModel%(ndim)id = PYB11TemplateClass(StrengthModel, template_parameters="%(Dimension)s")

LinearPolynomialEquationOfState%(ndim)id = PYB11TemplateClass(LinearPolynomialEquationOfState, template_parameters="%(Dimension)s")
GruneisenEquationOfState%(ndim)id = PYB11TemplateClass(GruneisenEquationOfState, template_parameters="%(Dimension)s")
OsborneEquationOfState%(ndim)id = PYB11TemplateClass(OsborneEquationOfState, template_parameters="%(Dimension)s")
TillotsonEquationOfState%(ndim)id = PYB11TemplateClass(TillotsonEquationOfState, template_parameters="%(Dimension)s")
MurnaghanEquationOfState%(ndim)id = PYB11TemplateClass(MurnaghanEquationOfState, template_parameters="%(Dimension)s")

SingularityEquationOfState%(ndim)id = PYB11TemplateClass(SingularityEquationOfState, template_parameters="Dim<%(ndim)i>")
SingularityGammaLawGas%(ndim)id = PYB11TemplateClass(SingularityGammaLawGas, template_parameters="Dim<%(ndim)i>")
SingularityExtendedVinet%(ndim)id = PYB11TemplateClass(SingularityExtendedVinet, template_parameters="Dim<%(ndim)i>")
SingularitySpiner%(ndim)id = PYB11TemplateClass(SingularitySpiner, template_parameters="Dim<%(ndim)i>")
SingularitySpinerT%(ndim)id = PYB11TemplateClass(SingularitySpinerT, template_parameters="Dim<%(ndim)i>")
SingularityEOSPAC%(ndim)id = PYB11TemplateClass(SingularityEOSPAC, template_parameters="Dim<%(ndim)i>")
    
NullStrength%(ndim)id = PYB11TemplateClass(NullStrength, template_parameters="%(Dimension)s")
ConstantStrength%(ndim)id = PYB11TemplateClass(ConstantStrength, template_parameters="%(Dimension)s")
SteinbergGuinanStrength%(ndim)id = PYB11TemplateClass(SteinbergGuinanStrength, template_parameters="%(Dimension)s")
JohnsonCookStrength%(ndim)id = PYB11TemplateClass(JohnsonCookStrength, template_parameters="%(Dimension)s")
CollinsStrength%(ndim)id = PYB11TemplateClass(CollinsStrength, template_parameters="%(Dimension)s")
iSALEROCKStrength%(ndim)id = PYB11TemplateClass(iSALEROCKStrength, template_parameters="%(Dimension)s")

PhysicsEvolvingMaterialLibrary%(ndim)id = PYB11TemplateClass(PhysicsEvolvingMaterialLibrary, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
