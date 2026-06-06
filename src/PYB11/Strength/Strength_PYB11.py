"""
Spheral Strength module.

Provides strength modeling helpers.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Strength/SolidFieldNames.hh"',
                  '"Geometry/Dimension.hh"',
                  '"Geometry/GeomPlane.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# SolidFieldNames
#-------------------------------------------------------------------------------
class SolidFieldNames:

    deviatoricStress = PYB11readonly(static=True, returnpolicy="copy")
    deviatoricStressTT = PYB11readonly(static=True, returnpolicy="copy")
    plasticStrain = PYB11readonly(static=True, returnpolicy="copy")
    plasticStrainRate = PYB11readonly(static=True, returnpolicy="copy")
    scalarDamage = PYB11readonly(static=True, returnpolicy="copy")
    tensorDamage = PYB11readonly(static=True, returnpolicy="copy")
    damageCoupling = PYB11readonly(static=True, returnpolicy="copy")
    strain = PYB11readonly(static=True, returnpolicy="copy")
    strainTensor = PYB11readonly(static=True, returnpolicy="copy")
    effectiveStrainTensor = PYB11readonly(static=True, returnpolicy="copy")
    bulkModulus = PYB11readonly(static=True, returnpolicy="copy")
    shearModulus = PYB11readonly(static=True, returnpolicy="copy")
    YoungsModulus = PYB11readonly(static=True, returnpolicy="copy")
    longitudinalSoundSpeed = PYB11readonly(static=True, returnpolicy="copy")
    yieldStrength = PYB11readonly(static=True, returnpolicy="copy")
    flaws = PYB11readonly(static=True, returnpolicy="copy")
    numFlaws = PYB11readonly(static=True, returnpolicy="copy")
    minFlaw = PYB11readonly(static=True, returnpolicy="copy")
    maxFlaw = PYB11readonly(static=True, returnpolicy="copy")
    initialVolume = PYB11readonly(static=True, returnpolicy="copy")
    randomGenerator = PYB11readonly(static=True, returnpolicy="copy")
    porositySolidDensity = PYB11readonly(static=True, returnpolicy="copy")
    porosityAlpha = PYB11readonly(static=True, returnpolicy="copy")
    porosityStrain = PYB11readonly(static=True, returnpolicy="copy")
    porosityAlpha0 = PYB11readonly(static=True, returnpolicy="copy")
    porosityc0 = PYB11readonly(static=True, returnpolicy="copy")
    fDSjutzi = PYB11readonly(static=True, returnpolicy="copy")
    fragmentIDs = PYB11readonly(static=True, returnpolicy="copy")
    particleTypes = PYB11readonly(static=True, returnpolicy="copy")
    meltSpecificEnergy = PYB11readonly(static=True, returnpolicy="copy")
    damagedPressure = PYB11readonly(static=True, returnpolicy="copy")
