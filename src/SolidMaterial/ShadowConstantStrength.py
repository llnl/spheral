#-------------------------------------------------------------------------------
# ShadowConstantStrength
#
# Provides convenient constructors for the ConstantStrength model using the canned
# values in MaterialPropertiesLib.py.
#-------------------------------------------------------------------------------
import types
from SpheralCompiledPackages import *
from MaterialPropertiesLib import SpheralMaterialPropertiesLib
from CaptureStdout import helpString
from spheralDimensions import spheralDimensions

dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Define a string providing the help for building a ConstantStrength.
#-------------------------------------------------------------------------------
expectedUsageString = """
ConstantStrength can be constructed one of two ways:

1.  Using canned material values stored in our internal data base.  Expected arguments:
        materialName        : Label for the material in data base

2.  You can directly set the strength parameters explicitly, as
        mu0                 : shear modulus
        Y0                  : yield strength for plastic yielding
"""

#-------------------------------------------------------------------------------
# The generic factory function, where you pass in the dimension specific 
# ConstantStrength constructor.
# This one is for internal use only -- people will actually call the dimension
# specific front-ends at the end of this script.
#-------------------------------------------------------------------------------
def _ConstantStrengthFactory(ndim):
    CXXConstantStrength = eval("ConstantStrength{}d".format(ndim))

    class ConstantStrength(CXXConstantStrength):

        def __init__(self,
                     *args, 
                     **kwargs):
            
            # The arguments that need to be passed to this method.
            expectedArgs = ["materialName", "units"]
            optionalKwArgs = {"mu0" : None,
                              "Y0" : None,
                              "muD" : 0.0,
                              "YD" : 0.0}

            # The base units for parameters in this file.
            cgs = PhysicalConstants(0.01,    # Length in m
                                    0.001,   # Mass in kg
                                    1.0)     # Time in sec

            # What sort of information did the user pass in?
            if ("materialName" in kwargs or 
                len(args) > 0 and type(args[0]) is str):

                # It looks like the user is trying to use one of the libarary canned values.
                # Evaluate the arguments to the method.
                if (len(args) > len(expectedArgs) or 
                    (len(args) + len(kwargs) < len(expectedArgs))): # insist on formal mandatory arguments 
                    raise ValueError(expectedUsageString)

                # Check for any invalid keywords
                for arg in kwargs: # deal with optional args
                    if arg not in (expectedArgs + list(optionalKwArgs.keys())):
                        raise ValueError(expectedUsageString)

                # Set the arguments dictionary
                dargs = {expectedArgs[i] : args[i] for i in range(len(args))}          # Mandatory args
                dargs.update(optionalKwArgs)
                dargs.update(kwargs)
                ARGS = types.SimpleNamespace(**dargs)

                # Check that the caller specified a valid material label.
                mat = ARGS.materialName.lower()
                if mat not in SpheralMaterialPropertiesLib:
                    raise ValueError("You must specify one of %s" % str(list(SpheralMaterialPropertiesLib.keys())))
                if ("mu0" not in SpheralMaterialPropertiesLib[mat] or
                    "Y0"  not in SpheralMaterialPropertiesLib[mat]):
                    raise ValueError("The material %s does not provide strength paramters." % materialName)

                # Extract the parameters for this material.
                if ARGS.mu0 is None:
                    ARGS.mu0 = SpheralMaterialPropertiesLib[mat]["mu0"]
                if ARGS.Y0 is None:
                    ARGS.Y0 = SpheralMaterialPropertiesLib[mat]["Y0"]
                if ARGS.muD is None:
                    ARGS.muD = ARGS.mu0
                if ARGS.YD is None:
                    ARGS.YD = ARGS.Y0

                # Figure out the conversions to the requested units.
                lconv = cgs.unitLengthMeters / ARGS.units.unitLengthMeters
                mconv = cgs.unitMassKg / ARGS.units.unitMassKg
                tconv = cgs.unitTimeSec / ARGS.units.unitTimeSec
                rhoConv = mconv/(lconv*lconv*lconv)
                Pconv = mconv/(lconv*tconv*tconv)
                specificEconv = (lconv/tconv)**2

                # Build the arguments for constructing the ConstantStrength.
                passargs = [ARGS.mu0 * Pconv,
                            ARGS.Y0 * Pconv]
                passkwargs = {}

            else:

                # Just pass through the arguments.
                passargs = args
                passkwargs = kwargs

            # Invoke the C++ constructor
            CXXConstantStrength.__init__(self, *tuple(passargs), **passkwargs)
            return

    return ConstantStrength

#-------------------------------------------------------------------------------
# Create the dimension specific ConstantStrength factories.  These are the ones
# you actually use.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec(f"ConstantStrength{ndim}d = _ConstantStrengthFactory({ndim})")
