from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

def FSISPH(dataBase,
           W,
           Q = None,
           slides=None,
           cfl = 0.35,
           surfaceForceCoefficient=0.0,
           densityStabilizationCoefficient=0.1, 
           specificThermalEnergyDiffusionCoefficient=0.1, 
           xsphCoefficient=0.0,
           interfaceMethod=HLLCInterface, 
           kernelAveragingMethod = NeverAverageKernels,      
           sumDensityNodeLists=[],
           useVelocityMagnitudeForDt = False,
           compatibleEnergyEvolution = True,
           evolveTotalEnergy = False,
           linearCorrectGradients = True,
           planeStrain = None,
           decoupleDamagedMaterial = True,
           interfacePmin = 0.0,
           interfaceNeighborAngleThreshold=0.707,
           HUpdate = IdealH,
           densityUpdate = FSISumMassDensity,
           epsTensile = 0.0,
           nTensile = 4.0,
           xmin = (-1e100, -1e100, -1e100),
           xmax = ( 1e100,  1e100,  1e100),
           etaMinAxis = 0.1,
           ASPH = False,
           smoothingScaleMethod = None):

    ######################################################################
    # some of these parameters are inactive and possible on there was out.
    # strengthInDamage and damageRelieveRubble are old switches and are not
    # implemented in the code.
    ######################################################################
    
    if compatibleEnergyEvolution and evolveTotalEnergy:
        raise RuntimeError("compatibleEnergyEvolution and evolveTotalEnergy are incompatible")

    if not planeStrain is None:
        print("FSISPH Warning: planeStrain option is deprecated and has no effect")

    # create the map nodelist --> index
    nodeLists = dataBase.nodeLists
    nodeListMap = {}
    for i in range(dataBase.numNodeLists):          
        nodeListMap[nodeLists[i]]=i

    # for now use int 1/0 to indicate sum density or not
    sumDensitySwitch = [0]*dataBase.numNodeLists
    for nodeList in sumDensityNodeLists:
        i = nodeListMap[nodeList]
        sumDensitySwitch[i]=1
    sumDensitySwitch = vector_of_int(sumDensitySwitch)
            
    # We use the provided DataBase to sniff out what sort of NodeLists are being
    # used, and based on this determine which SPH object to build.
    ndim = dataBase.nDim
    nfluid = dataBase.numFluidNodeLists
    nsolid = dataBase.numSolidNodeLists
    if nsolid > 0 and nsolid != nfluid:
        print("FSISPH Error: you have provided both solid and fluid NodeLists, which is currently not supported.")
        print("              If you want some fluids active, provide SolidNodeList without a strength option specfied,")
        print("              which will result in fluid behaviour for those nodes.")
        raise RuntimeError("Cannot mix solid and fluid NodeLists.")

    # Decide on the hydro object.
    if GeometryRegistrar.coords() == CoordinateType.Spherical:
        assert ndim == 1
        raise RuntimeError("not implemented for spherical coordinates")
    if GeometryRegistrar.coords() == CoordinateType.RZ:
        assert ndim == 2
        if nsolid > 0:
            Constructor = SolidFSISPHRZ
        else:
            raise RuntimeError("currently only implemented for solid nodelists")
    else:
        # Cartesian ---------------------------------
        if nsolid > 0:
            Constructor = eval("SolidFSISPH%id" % ndim)
        else:
            raise RuntimeError("currently only implemented for solid nodelists")


    # Artificial viscosity.
    if not Q:
        Cl = 2.0*(dataBase.maxKernelExtent/2.0)
        Cq = 8.0*(dataBase.maxKernelExtent/2.0)**2
        Q = eval("LimitedMonaghanGingoldViscosity%id(Clinear=%g, Cquadratic=%g, kernel=W)" % (ndim, Cl, Cq))

    # slide surfaces.
    if not slides:
        contactTypes = vector_of_int([0]*(dataBase.numNodeLists**2))
        slides = eval("SlideSurface%id(dataBase,contactTypes)" % ndim)

    # Build the constructor arguments
    xmin = (ndim,) + xmin
    xmax = (ndim,) + xmax
    kwargs = {"dataBase" : dataBase,
              "Q" : Q,
              "slides" : slides,
              "W" : W,
              "cfl" : cfl,
              "surfaceForceCoefficient" : surfaceForceCoefficient,
              "densityStabilizationCoefficient" : densityStabilizationCoefficient,      
              "specificThermalEnergyDiffusionCoefficient" : specificThermalEnergyDiffusionCoefficient,        
              "xsphCoefficient" : xsphCoefficient,
              "interfaceMethod" : interfaceMethod,
              "kernelAveragingMethod" : kernelAveragingMethod, 
              "sumDensityNodeLists" : sumDensitySwitch,
              "useVelocityMagnitudeForDt" : useVelocityMagnitudeForDt,
              "compatibleEnergyEvolution" : compatibleEnergyEvolution,
              "evolveTotalEnergy" : evolveTotalEnergy,
              "linearCorrectGradients" : linearCorrectGradients,
              "decoupleDamagedMaterial" : decoupleDamagedMaterial,
              "interfacePmin" : interfacePmin,
              "interfaceNeighborAngleThreshold" : interfaceNeighborAngleThreshold,
              "densityUpdate" : densityUpdate,
              "epsTensile" : epsTensile,
              "nTensile" : nTensile,
              "xmin" : eval("Vector%id(%g, %g, %g)" % xmin),
              "xmax" : eval("Vector%id(%g, %g, %g)" % xmax)}

    # Build and return the thing.
    result = Constructor(**kwargs)
    result.Q = Q
    result.slides = slides
    
    # Add the Q as a sub-package (to run before the hydro)
    result.prependSubPackage(Q)

    # Smoothing scale update
    if smoothingScaleMethod is None:
        if ASPH:
            if isinstance(ASPH, str) and ASPH.upper() == "CLASSIC":
                smoothingScaleMethod = eval(f"ASPHClassicSmoothingScale{ndim}d({HUpdate}, W)")
            else:
                smoothingScaleMethod = eval(f"ASPHSmoothingScale{ndim}d({HUpdate}, W)")
        else:
            smoothingScaleMethod = eval(f"SPHSmoothingScale{ndim}d({HUpdate}, W)")
    result._smoothingScaleMethod = smoothingScaleMethod
    result.appendSubPackage(smoothingScaleMethod)

    if GeometryRegistrar.coords() == CoordinateType.RZ:
        result.zaxisBC = AxisBoundaryRZ(etaMinAxis)
        result.appendBoundary(result.zaxisBC)

    return result

#-------------------------------------------------------------------------------
# Provide shorthand names for AFSISPH
#-------------------------------------------------------------------------------
def AFSISPH(*args, **kwargs):
    kwargs.update({"ASPH" : True})
    return SPH(*args, **kwargs)
