import sys
from SpheralCompiledPackages import *
from MaterialPropertiesLib import SpheralMaterialPropertiesLib

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# ProbabilisticDamageModel factory function
#-------------------------------------------------------------------------------
def PDMFactory(ndim):
    IntField = eval("IntField{}d".format(ndim))
    CXXProbabilisticDamageModel = eval("ProbabilisticDamageModel{}d".format(ndim))
    class ProbabilisticDamageModel(CXXProbabilisticDamageModel):

        def __init__(self, *args, **kwargs):
            """
ProbabilisticDamageModel is constructed with the following arguments (any 
default values listed in parens):

        materialName            : (optional) label for the material in data base to
                                  lookup the Weibull parameters
        units                   : (optional) scale kWeibull/mWeibull lookup from
                                  material name to given units
        nodeList                : (required) the FluidNodeList this model should be 
                                  applied to
        kernel                  : (required) the interpolation kernel to use
        kWeibull                : the "k" Weibull constant -- can be looked up from 
                                  materialName or provided
        mWeibull                : the "m" Weibull constant -- can be looked up from 
                                  materialName or provided
        seed                    : (48927592) random number seed for flaw generation
        minFlawsPerNode         : (1) the minimum number of flaws to seed on a point
        crackGrowthMultiplier   : (0.4) crack growth rate in units of longitudinal
                                  sound speed
        volumeMultiplier        : (1.0) Multiplies per node volume, useful for 
                                  reduced dimensional (1D or 2D) problems.
        damageCouplingAlgorithm : (PairMaxDamage) how should damaged points couple
        strainAlgorithm         : (PseudoPlasticStrain) how to define the strain
        damageInCompression     : (False) optionally damage in compression as well 
                                  as tension
        criticalDamageThreshold : (4.0) prevent any nodes where Trace(D_i) exceeds 
                                  criticalDamageThreshold from setting timestep
        mask                    : (1 on all points) a field of flags: a node with 
                                  zero implies no flaws (and therefore no damage) 
                                  on that point
"""

            # The material library values are in CGS, so build a units object for 
            # conversions.
            cgs = PhysicalConstants(0.01,   # unit length in m
                                    0.001,  # unit mass in kg
                                    1.0)    # unit time in sec

            # Arguments needed to build the damage model.
            damage_kwargs = {"nodeList"                 : None,
                             "kernel"                   : None,
                             "kWeibull"                 : None,
                             "mWeibull"                 : None,
                             "seed"                     : 48927592,
                             "minFlawsPerNode"          : 1,
                             "crackGrowthMultiplier"    : 0.4,
                             "volumeMultiplier"         : 1.0,
                             "damageCouplingAlgorithm"  : PairMaxDamage,
                             "strainAlgorithm"          : PseudoPlasticStrain,
                             "damageInCompression"      : False,
                             "criticalDamageThreshold"  : 4.0,
                             "mask"                     : None}

            # Extra arguments for our convenient constructor.
            convenient_kwargs = {"materialName"          : None,
                                 "units"                 : None}

            # Check the input arguments.
            validKeys = list(damage_kwargs.keys()) + list(convenient_kwargs.keys())
            for argname in kwargs:
                if not argname in validKeys:
                    raise ValueError("ERROR: argument {} not a valid option.\n".format(argname))

            # Did the user try any convenient constructor operations?
            if ((len(args) > 0 and type(args[0]) == str) or
                "materialName" in kwargs):
                if len(args) > 0 and type(args[0]) == str:
                    materialName = args[0]
                    args = args[1:]
                else:
                    materialName = kwargs["materialName"]
                    del kwargs["materialName"]
                if not materialName in SpheralMaterialPropertiesLib:
                    raise ValueError("ERROR: material {} is not in the library of material values.\n".format(materialName))
                matprops = SpheralMaterialPropertiesLib[materialName]
                if not ("kWeibull" in matprops and "mWeibull" in matprops):
                    raise ValueError("ERROR : material {} does not provide the required values for kWeibull and mWeibull.\n".format(materialName))
                damage_kwargs["kWeibull"] = matprops["kWeibull"]
                damage_kwargs["mWeibull"] = matprops["mWeibull"]

                # Any attempt to specify units?
                if "units" in kwargs:
                    units = kwargs["units"]
                    damage_kwargs["kWeibull"] *= (cgs.unitLengthMeters/units.unitLengthMeters)**3
                    del kwargs["units"]
                elif len(args) > 1 and isinstance(args[1], PhysicalConstants):
                    units = args[1]
                    damage_kwargs["kWeibull"] *= (cgs.unitLengthMeters/units.unitLengthMeters)**3
                    args = args[1:]

            # Process remaining user arguments.
            kwarg_order = ["nodeList", "kernel", "kWeibull", "mWeibull", "seed", "minFlawsPerNode", "crackGrowthMultiplier",
                           "volumeMultiplier", "damageCouplingAlgorithm", "strainAlgorithm", "damageInCompression",
                           "criticalDamageThreshold", "mask"]
            for iarg, argval in enumerate(args):
                damage_kwargs[kwarg_order[iarg]] = argval

            # Process any keyword arguments.  Note we already removed any deprecated keywords.
            for argname in kwargs:
                if argname in damage_kwargs:
                    damage_kwargs[argname] = kwargs[argname]
                else:
                    raise ValueError("ERROR : unknown kwarg {}.\n".format(argname))

            # If no mask was provided, deafult to all points active
            if damage_kwargs["mask"] is None:
                damage_kwargs["mask"] = IntField("damage mask", damage_kwargs["nodeList"], 1)

            # # If no minFlawsPerNode was specified, set it to a fraction of log(N_points)
            # if damage_kwargs["minFlawsPerNode"] is None:
            #     damage_kwargs["minFlawsPerNode"] = max(1, int(log(mpi.allreduce(damage_kwargs["nodeList"].numInternalNodes, mpi.SUM))))

            # Invoke the parent constructor
            CXXProbabilisticDamageModel.__init__(self, **damage_kwargs)

    return ProbabilisticDamageModel

#-------------------------------------------------------------------------------
# Create the different dimension implementations.
#-------------------------------------------------------------------------------
for ndim in dims:
    exec("ProbabilisticDamageModel{ndim}d = PDMFactory({ndim})".format(ndim=ndim))
