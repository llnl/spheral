from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# linear spring helper class
#-------------------------------------------------------------------------------
def LinearSpringDEM(dataBase,
                    normalSpringConstant,
                    normalRestitutionCoefficient,
                    tangentialSpringConstant,
                    tangentialRestitutionCoefficient,
                    dynamicFrictionCoefficient,
                    staticFrictionCoefficient,
                    rollingFrictionCoefficient,
                    torsionalFrictionCoefficient,
                    cohesiveTensileStrength = 0.0,
                    shapeFactor = 0.0,
                    stepsPerCollision = 25,
                    enableFastTimeStepping = False,
                    xmin = (-1e100, -1e100, -1e100),
                    xmax = ( 1e100,  1e100,  1e100),
                    normalRestitutionCoefficientParticleBoundary = None,
                    tangentialRestitutionCoefficientParticleBoundary = None,
                    dynamicFrictionCoefficientParticleBoundary = None,
                    staticFrictionCoefficientParticleBoundary = None,
                    rollingFrictionCoefficientParticleBoundary = None,
                    torsionalFrictionCoefficientParticleBoundary = None):

    assert dataBase.numDEMNodeLists == dataBase.numNodeLists, "all nodelists must be dem nodelists"
    assert stepsPerCollision > 1, "stepsPerCollision too low, recommended is 25-50"
    assert cohesiveTensileStrength >= 0, "cohesiveTensileStrength must be positive"
    assert normalSpringConstant >= 0, "normalSpringConstant must be positive"
    assert normalRestitutionCoefficient >= 0 and normalRestitutionCoefficient <= 1, "normalSpringConstant must be between 1 and 0"
    assert tangentialSpringConstant >= 0, "normalSpringConstant must be positive"
    assert tangentialRestitutionCoefficient >= 0 and tangentialRestitutionCoefficient <= 1, "normalSpringConstant must be between 1 and 0"
    assert staticFrictionCoefficient >= 0, "staticFrictionCoefficient must be positive"
    assert dynamicFrictionCoefficient >= 0, "dynamicFrictionCoefficient must be positive"
    assert rollingFrictionCoefficient >= 0, "rollingFrictionCoefficient must be positive"
    assert torsionalFrictionCoefficient >= 0, "torsionalFrictionCoefficient must be positive"
    assert isinstance(enableFastTimeStepping,bool)
    
    if normalRestitutionCoefficientParticleBoundary is None:
        normalRestitutionCoefficientParticleBoundary = normalRestitutionCoefficient
    if tangentialRestitutionCoefficientParticleBoundary is None:
        tangentialRestitutionCoefficientParticleBoundary = tangentialRestitutionCoefficient
    if dynamicFrictionCoefficientParticleBoundary is None:
        dynamicFrictionCoefficientParticleBoundary = dynamicFrictionCoefficient
    if staticFrictionCoefficientParticleBoundary is None:
        staticFrictionCoefficientParticleBoundary = staticFrictionCoefficient
    if rollingFrictionCoefficientParticleBoundary is None:
        rollingFrictionCoefficientParticleBoundary = rollingFrictionCoefficient
    if torsionalFrictionCoefficientParticleBoundary is None:
        torsionalFrictionCoefficientParticleBoundary = torsionalFrictionCoefficient

    assert normalRestitutionCoefficientParticleBoundary >= 0 and normalRestitutionCoefficientParticleBoundary <= 1, "normalRestitutionCoefficientParticleBoundary must be between 1 and 0"
    assert tangentialRestitutionCoefficientParticleBoundary >= 0 and tangentialRestitutionCoefficientParticleBoundary <= 1, "tangentialRestitutionCoefficientParticleBoundary must be between 1 and 0"
    assert dynamicFrictionCoefficientParticleBoundary >= 0, "dynamicFrictionCoefficientParticleBoundary must be positive"
    assert staticFrictionCoefficientParticleBoundary >= 0, "staticFrictionCoefficientParticleBoundary must be positive"
    assert rollingFrictionCoefficientParticleBoundary >= 0, "rollingFrictionCoefficientParticleBoundary must be positive"
    assert torsionalFrictionCoefficientParticleBoundary >= 0, "torsionalFrictionCoefficientParticleBoundary must be positive"

    #if (stepsPerCollision < 10) print("WARNING: stepsPerCollision is very low, recommended is 25-50")

    # we might want to allow the user to set less parameters with reasonable defaults
    #if tangentialSpringConstant is None:
    #    tangentialSpringConstant = normalSpringConstant * 2.0/7.0
    #if tangentialRestitutionCoefficient is None:
    #    tangentialRestitutionCoefficient = normalRestitutionCoefficient

    ndim = dataBase.nDim

    Constructor = eval("LinearSpringDEM%id" % ndim)

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"dataBase" : dataBase,
              "normalSpringConstant" : normalSpringConstant,
              "normalRestitutionCoefficient" : normalRestitutionCoefficient,
              "tangentialSpringConstant" : tangentialSpringConstant,
              "tangentialRestitutionCoefficient" : tangentialRestitutionCoefficient,
              "dynamicFrictionCoefficient" : dynamicFrictionCoefficient,
              "staticFrictionCoefficient" : staticFrictionCoefficient,
              "rollingFrictionCoefficient" : rollingFrictionCoefficient,
              "torsionalFrictionCoefficient" : torsionalFrictionCoefficient,
              "normalRestitutionCoefficientParticleBoundary" : normalRestitutionCoefficientParticleBoundary,
              "tangentialRestitutionCoefficientParticleBoundary" : tangentialRestitutionCoefficientParticleBoundary,
              "dynamicFrictionCoefficientParticleBoundary" : dynamicFrictionCoefficientParticleBoundary,
              "staticFrictionCoefficientParticleBoundary" : staticFrictionCoefficientParticleBoundary,
              "rollingFrictionCoefficientParticleBoundary" : rollingFrictionCoefficientParticleBoundary,
              "torsionalFrictionCoefficientParticleBoundary" : torsionalFrictionCoefficientParticleBoundary,
              "cohesiveTensileStrength" : cohesiveTensileStrength,
              "shapeFactor" : shapeFactor,
              "stepsPerCollision" : stepsPerCollision,
              "enableFastTimeStepping" : enableFastTimeStepping,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    # Build and return the thing.
    result = Constructor(**kwargs)

    return result

#-------------------------------------------------------------------------------
# convienence function that defaults to Linear Spring DEM
#-------------------------------------------------------------------------------
def DEM(dataBase,
        normalSpringConstant,
        normalRestitutionCoefficient,
        tangentialSpringConstant,
        tangentialRestitutionCoefficient,
        dynamicFrictionCoefficient,
        staticFrictionCoefficient,
        rollingFrictionCoefficient,
        torsionalFrictionCoefficient,
        cohesiveTensileStrength=0.0,
        shapeFactor=0.0,
        stepsPerCollision = 25,
        enableFastTimeStepping = True,
        xmin = (-1e100, -1e100, -1e100),
        xmax = ( 1e100,  1e100,  1e100),
        normalRestitutionCoefficientParticleBoundary = None,
        tangentialRestitutionCoefficientParticleBoundary = None,
        dynamicFrictionCoefficientParticleBoundary = None,
        staticFrictionCoefficientParticleBoundary = None,
        rollingFrictionCoefficientParticleBoundary = None,
        torsionalFrictionCoefficientParticleBoundary = None):
    return LinearSpringDEM(dataBase = dataBase,
                           normalSpringConstant = normalSpringConstant,
                           normalRestitutionCoefficient = normalRestitutionCoefficient,
                           tangentialSpringConstant = tangentialSpringConstant,
                           tangentialRestitutionCoefficient = tangentialRestitutionCoefficient,
                           dynamicFrictionCoefficient = dynamicFrictionCoefficient,
                           staticFrictionCoefficient = staticFrictionCoefficient,
                           rollingFrictionCoefficient = rollingFrictionCoefficient,
                           torsionalFrictionCoefficient = torsionalFrictionCoefficient,
                           cohesiveTensileStrength = cohesiveTensileStrength,
                           shapeFactor = shapeFactor,
                           stepsPerCollision = stepsPerCollision,
                           enableFastTimeStepping = enableFastTimeStepping,
                           xmin = xmin,
                           xmax = xmax,
                           normalRestitutionCoefficientParticleBoundary = normalRestitutionCoefficientParticleBoundary,
                           tangentialRestitutionCoefficientParticleBoundary = tangentialRestitutionCoefficientParticleBoundary,
                           dynamicFrictionCoefficientParticleBoundary = dynamicFrictionCoefficientParticleBoundary,
                           staticFrictionCoefficientParticleBoundary = staticFrictionCoefficientParticleBoundary,
                           rollingFrictionCoefficientParticleBoundary = rollingFrictionCoefficientParticleBoundary,
                           torsionalFrictionCoefficientParticleBoundary = torsionalFrictionCoefficientParticleBoundary)
