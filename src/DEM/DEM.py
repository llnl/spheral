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
                    xmax = ( 1e100,  1e100,  1e100)):

    assert dataBase.numDEMNodeLists == dataBase.numNodeLists, "all nodelists must be dem nodelists"
    assert stepsPerCollision > 1, "stepsPerCollision too low, recommended is 25-50"
    assert cohesiveTensileStrength >= 0, "cohesiveTensileStrength must be positive"
    assert normalSpringConstant >= 0, "normalSpringConstant must be positive"
    assert normalRestitutionCoefficient >= 0 and normalRestitutionCoefficient <= 1, "normalSpringConstant must be between 1 and 0"
    assert tangentialSpringConstant >= 0, "normalSpringConstant must be positive"
    assert tangentialRestitutionCoefficient >= 0 and tangentialRestitutionCoefficient <= 1, "normalSpringConstant must be between 1 and 0"
    assert dynamicFrictionCoefficient >= 0, "dynamicFrictionCoefficient must be positive"
    assert dynamicFrictionCoefficient >= 0, "dynamicFrictionCoefficient must be positive"
    assert rollingFrictionCoefficient >= 0, "rollingFrictionCoefficient must be positive"
    assert torsionalFrictionCoefficient >= 0, "torsionalFrictionCoefficient must be positive"
    assert isinstance(enableFastTimeStepping,bool)

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
        xmax = ( 1e100,  1e100,  1e100)):
    return LinearSpringDEM(dataBase,
                           normalSpringConstant,
                           normalRestitutionCoefficient,
                           tangentialSpringConstant,
                           tangentialRestitutionCoefficient,
                           dynamicFrictionCoefficient,
                           staticFrictionCoefficient,
                           rollingFrictionCoefficient,
                           torsionalFrictionCoefficient,
                           cohesiveTensileStrength,
                           shapeFactor,
                           stepsPerCollision,
                           enableFastTimeStepping,
                           xmin,
                           xmax)
