#-------------------------------------------------------------------------------
# UnitConversion
#-------------------------------------------------------------------------------
from PYB11Generator import *

class UnitConversion:
    """
    Convert between unit systems
    """

    #...........................................................................
    # Constructor
    def pyinit(self,
               unitsFrom = "const PhysicalConstants&",
               unitsTo = "const PhysicalConstants&"):
        "Convert from unitsFrom to unitsTo"
        return

    #...........................................................................
    # Properties
    length = PYB11property("double", "length")
    mass = PYB11property("double", "mass")
    time = PYB11property("double", "time")
    temperature = PYB11property("double", "temperature")
    charge = PYB11property("double", "charge")
    massDensity = PYB11property("double", "massDensity")
    specificEnergy = PYB11property("double", "specificEnergy")
    energyDensity = PYB11property("double", "energyDensity")
    pressure = PYB11property("double", "pressure")
    bulkModulus = PYB11property("double", "bulkModulus")
    specificHeat = PYB11property("double", "specificHeat")
    entropy = PYB11property("double", "entropy")
    speed = PYB11property("double", "speed")
    volume = PYB11property("double", "volume")
