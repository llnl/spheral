from PYB11Generator import *

@PYB11holder("std::shared_ptr")
class HypreOptions:
    def pyinit(self):
        "Holds the options for Hypre preconditioners and solver"

    quitIfDiverged = PYB11readwrite()
    
    logLevel = PYB11readwrite()
    printLevel = PYB11readwrite()

    sumDuplicates = PYB11readwrite()
    addToValues = PYB11readwrite()
    
    saveLinearSystem = PYB11readwrite()
    printIterations = PYB11readwrite()
    maxNumberOfIterations = PYB11readwrite()
    toleranceL2 = PYB11readwrite()
    useRobustTolerance = PYB11readwrite()
    absoluteTolerance = PYB11readwrite()

    kDim = PYB11readwrite()
    minIters = PYB11readwrite()
    
    HyprePreconditionerType = PYB11enum(("NoPreconditioner",
                                         "AMGPreconditioner",
                                         "ILUPreconditioner"),
                                         export_values = True,
                                         doc = "Preconditioners available for Hypre")
    preconditionerType = PYB11readwrite()
    
    measure_type = PYB11readwrite()
    pcTol = PYB11readwrite()
    minItersAMG = PYB11readwrite()
    maxItersAMG = PYB11readwrite()
    cycleType = PYB11readwrite()
    coarsenTypeAMG = PYB11readwrite()
    strongThresholdAMG = PYB11readwrite()
    maxRowSumAMG = PYB11readwrite()
    interpTypeAMG = PYB11readwrite()
    aggNumLevelsAMG = PYB11readwrite()
    aggInterpTypeAMG = PYB11readwrite()
    pMaxElmtsAMG = PYB11readwrite()
    truncFactorAMG = PYB11readwrite()
    printLevelAMG = PYB11readwrite()
    logLevelAMG = PYB11readwrite()
    relaxWeightAMG = PYB11readwrite()
    relaxTypeAMG = PYB11readwrite()
    relaxTypeCoarseAMG = PYB11readwrite()
    cycleNumSweepsAMG = PYB11readwrite()
    cycleNumSweepsCoarseAMG = PYB11readwrite()
    maxLevelsAMG = PYB11readwrite()

    useILUT = PYB11readwrite()
    factorLevelILU = PYB11readwrite()
    rowScaleILU = PYB11readwrite()
    printLevelILU = PYB11readwrite()
    dropToleranceILU = PYB11readwrite()
    
    zeroOutNegativities = PYB11readwrite()

