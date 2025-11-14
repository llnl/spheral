from PYB11Generator import *
from SingularityEquationOfState import *

@PYB11template("Dimension")
class SingularityGammaLawGas(SingularityEquationOfState):
    def pyinit(self,
               gamma = "const double",
               mu = "const double",
               constants = "const PhysicalConstants&",
               level = ("const SingularityLimitLevel", "SingularityLimitLevel::NONE"),
               referenceDensity = ("const double", "1.0"),
               etamin = ("const double", "0.0"),
               etamax = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double",  "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor"),
               externalPressure = ("const double", "0.0")):
        "Gamma law gas"

@PYB11template("Dimension")
class SingularityExtendedVinet(SingularityEquationOfState):
    def pyinit(self,
               rho0 = "const double",
               T0 = "const double",
               B0 = "const double",
               BP0 = "const double",
               A0 = "const double",
               Cv0 = "const double",
               E0 = "const double",
               S0 = "const double",
               coeffs = "const std::vector<double>&",
               constants = "const PhysicalConstants&",
               level = ("const SingularityLimitLevel", "SingularityLimitLevel::NONE"),
               referenceDensity = ("const double", "1.0"),
               etamin = ("const double", "0.0"),
               etamax = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double",  "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor"),
               externalPressure = ("const double", "0.0")):
        "Extended Vinet"

# @PYB11template("Dimension")
# class SingularityMieGruneisenLinear(SingularityEquationOfState):
#     def pyinit(self,
#                rho0 = "const double",
#                T0 = "const double",
#                Cs = "const double",
#                s = "const double",
#                G0 = "const double",
#                Cv0 = "const double",
#                E0 = "const double",
#                S0 = "const double",
#                constants = "const PhysicalConstants&",
#                level = ("const SingularityLimitLevel", "SingularityLimitLevel::NONE"),
#                minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
#                maximumPressure = ("const double",  "std::numeric_limits<double>::max()"),
#                minimumPressureDamage = ("const double", "0.0"),
#                minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor"),
#                externalPressure = ("const double", "0.0")):
#         "Mie Gruneisen"

@PYB11template("Dimension")
class SingularitySpiner(SingularityEquationOfState):
    def pyinit(self,
               filename = "const std::string&",
               matid = "int",
               constants = "const PhysicalConstants&",
               level = ("const SingularityLimitLevel", "SingularityLimitLevel::CS"),
               referenceDensity = ("const double", "1.0"),
               etamin = ("const double", "0.0"),
               etamax = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double",  "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor"),
               externalPressure = ("const double", "0.0")):
        "Spiner"

@PYB11template("Dimension")
class SingularitySpinerT(SingularityEquationOfState):
    def pyinit(self,
               filename = "const std::string&",
               matid = "int",
               constants = "const PhysicalConstants&",
               level = ("const SingularityLimitLevel", "SingularityLimitLevel::CS"),
               referenceDensity = ("const double", "1.0"),
               etamin = ("const double", "0.0"),
               etamax = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double",  "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor"),
               externalPressure = ("const double", "0.0")):
        "Spiner with lookup in T"

@PYB11template("Dimension")
class SingularityEOSPAC(SingularityEquationOfState):
    def pyinit(self,
               matid = "int",
               constants = "const PhysicalConstants&",
               level = ("const SingularityLimitLevel", "SingularityLimitLevel::NONE"),
               referenceDensity = ("const double", "1.0"),
               etamin = ("const double", "0.0"),
               etamax = ("const double", "std::numeric_limits<double>::max()"),
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double",  "std::numeric_limits<double>::max()"),
               minimumPressureDamage = ("const double", "0.0"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor"),
               externalPressure = ("const double", "0.0")):
        "EOSPAC"

