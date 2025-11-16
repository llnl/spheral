from math import *

from NodeGeneratorBase import *

from Spheral import Vector1d, Vector2d, Vector3d, SymTensor1d, SymTensor2d, SymTensor3d
from SpheralTestUtilities import fuzzyEqual

#-------------------------------------------------------------------------------
# Wrapper Generator for DEM based on SPH generators
#-------------------------------------------------------------------------------
class GenerateDEMfromSPHGenerator1d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 W,
                 SPHGenerator,
                 DEMParticleGenerator=None,
                 particleRadius=None): 

        # set kernel equal extent to 2x particle radius
        # this is a dummy val for now. DEM package 
        # will re-adjust to appropriate H later.
        hOverR = 2.0/W.kernelExtent

        def constantRadiusFunc(position):
                return particleRadius

        # set up our initial radius
        #--------------------------------------------------
        if isinstance(particleRadius, list):
            radii = list(particleRadius)
            counter = {"i" : 0}
            def radiusFunc(position):
                r = radii[counter["i"]]
                counter["i"] += 1
                return r
            self.particleRadiusFunc = radiusFunc
        elif type(particleRadius) in [float,int]:
            self.particleRadiusFunc = constantRadiusFunc
        elif hasattr(particleRadius,"__call__"):
            self.particleRadiusFunc = particleRadius
        else:
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("WARNING: no particle radius specified, attempting to put something together...")
            print("         calculating min h from SPH generator...")
            hmin = 1e100
            for Hi in SPHGenerator.H:
                hi = Hi.Inverse().Trace() / SymTensor1d.nDimensions
                hmin = min(hmin,hi)
            print("         calculating node spacing from min h...")
            nodeSpacing = hmin/SPHGenerator.nNodePerh
            print("         calculating constant particle radius...")
            particleRadius = 0.48*nodeSpacing
            self.particleRadiusFunc = constantRadiusFunc
            print("         it worked, check your particle radii to make sure its what you wanted.")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            

        # create composite particle distribution from base
        #--------------------------------------------------
        if DEMParticleGenerator:
            self.x = []
            self.m = []
            self.particleRadius = []
            self.compositeParticleIndex = []
            for i in range(SPHGenerator.localNumNodes()):

                xi,mi,Ri = DEMParticleGenerator(SPHGenerator.x[i],
                                                SPHGenerator.H[i],
                                                SPHGenerator.m[i],
                                                self.particleRadiusFunc(Vector1d(SPHGenerator.x[i])))

                pIDi = [SPHGenerator.globalIDs[i]+1 for entry in Ri]
                self.x.extend(xi)
                self.m.extend(mi)
                self.particleRadius.extend(Ri)
                self.compositeParticleIndex.extend(pIDi)
                
        else:
            self.x = SPHGenerator.x
            self.m = SPHGenerator.m
            self.particleRadius = [self.particleRadiusFunc(Vector1d(SPHGenerator.x[i])) for i in range(len(self.x))]
            self.compositeParticleIndex = SPHGenerator.globalIDs
    
        self.H = [SymTensor1d(1.0/(Rj*hOverR)) for Rj in self.particleRadius]

        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.m, self.H, self.particleRadius, self.compositeParticleIndex)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localParticleRadius(self, i):
        return self.particleRadius[i]

    def globalParticleRadius(self, i):
        return self._globalValue(i, self.localParticleRadius)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localCompositeParticleIndex(self, i):
        return self.compositeParticleIndex[i]

    def globalCompositeParticleIndex(self, i):
        return self._globalValue(i, self.localCompositeParticleIndex)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        return Vector1d(self.x[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        return self.H[i]
    

#-------------------------------------------------------------------------------
# Wrapper Generator for DEM based on SPH generators
#-------------------------------------------------------------------------------
class GenerateDEMfromSPHGenerator2d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 W,
                 SPHGenerator,
                 DEMParticleGenerator=None,
                 particleRadius=None): 

        # set kernel equal extent to 2x particle radius
        # this is a dummy val for now. DEM package 
        # will re-adjust to appropriate H later.
        hOverR = 2.0/W.kernelExtent

        def constantRadiusFunc(position):
                return particleRadius

        # set up our initial radius
        #--------------------------------------------------
        if type(particleRadius) in [float,int]:
            self.particleRadiusFunc = constantRadiusFunc
        elif hasattr(particleRadius,"__call__"):
            self.particleRadiusFunc = particleRadius
        else:
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("WARNING: no particle radius specified, attempting to put something together...")
            print("         calculating min h from SPH generator...")
            hmin = 1e100
            for Hi in SPHGenerator.H:
                hi = Hi.Inverse().Trace() / SymTensor2d.nDimensions
                hmin = min(hmin,hi)
            print("         calculating node spacing from min h...")
            nodeSpacing = hmin/SPHGenerator.nNodePerh
            print("         calculating constant particle radius...")
            particleRadius = 0.48*nodeSpacing
            self.particleRadiusFunc = constantRadiusFunc
            print("         it worked, check your particle radii to make sure its what you wanted.")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            

        # create composite particle distribution from base
        #--------------------------------------------------
        if DEMParticleGenerator:
            self.x = []
            self.y = []
            self.m = []
            self.particleRadius = []
            self.compositeParticleIndex = []
            for i in range(SPHGenerator.localNumNodes()):

                xi,yi,mi,Ri = DEMParticleGenerator(SPHGenerator.x[i],
                                                   SPHGenerator.y[i],
                                                   SPHGenerator.H[i],
                                                   SPHGenerator.m[i],
                                                   self.particleRadiusFunc(Vector2d(SPHGenerator.x[i],SPHGenerator.y[i])))

                pIDi = [SPHGenerator.globalIDs[i]+1 for entry in Ri]
                self.x.extend(xi)
                self.y.extend(yi)
                self.m.extend(mi)
                self.particleRadius.extend(Ri)
                self.compositeParticleIndex.extend(pIDi)
                
        else:
            self.x = SPHGenerator.x
            self.y = SPHGenerator.y
            self.m = SPHGenerator.m
            self.particleRadius = [self.particleRadiusFunc(Vector2d(SPHGenerator.x[i],SPHGenerator.y[i])) for i in range(len(self.x))]
            self.compositeParticleIndex = SPHGenerator.globalIDs
    
        self.H = [SymTensor2d(1.0/(Rj*hOverR), 0.0, 0.0, 1.0/(Rj*hOverR)) for Rj in self.particleRadius]

        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.m, self.H, self.particleRadius, self.compositeParticleIndex)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localParticleRadius(self, i):
        return self.particleRadius[i]

    def globalParticleRadius(self, i):
        return self._globalValue(i, self.localParticleRadius)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localCompositeParticleIndex(self, i):
        return self.compositeParticleIndex[i]

    def globalCompositeParticleIndex(self, i):
        return self._globalValue(i, self.localCompositeParticleIndex)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        return Vector2d(self.x[i],self.y[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        return self.H[i]
    

#-------------------------------------------------------------------------------
# Wrapper Generator for DEM based on SPH generators
#-------------------------------------------------------------------------------
class GenerateDEMfromSPHGenerator3d(NodeGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 W,
                 SPHGenerator,
                 DEMParticleGenerator=None,
                 particleRadius=None): 

        # set kernel equal extent to 2x particle radius
        # this is a dummy val for now. DEM package 
        # will re-adjust to appropriate H later.
        hOverR = 2.0/W.kernelExtent

        def constantRadiusFunc(position):
            return particleRadius

        # set up our initial radius
        #--------------------------------------------------
        if type(particleRadius) in [float,int]:
            def constantRadiusFunc(position):
                return particleRadius
            self.particleRadiusFunc = constantRadiusFunc
        elif hasattr(particleRadius,"__call__"):
            self.particleRadiusFunc = constantRadiusFunc
        else:
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print("WARNING: no particle radius specified, attempting to put something together...")
            print("         calculating min h from SPH generator...")
            hmin = 1e100
            for Hi in SPHGenerator.H:
                hi = Hi.Inverse().Trace() / SymTensor3d.nDimensions
                hmin = min(hmin,hi)
            print("         calculating node spacing from min h...")
            nodeSpacing = hmin/SPHGenerator.nNodePerh
            print("         calculating constant particle radius...")
            particleRadius = 0.48*nodeSpacing
            self.particleRadiusFunc = constantRadiusFunc
            print("         it worked, check your particle radii to make sure its what you wanted.")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            

        # create composite particle distribution from base
        #--------------------------------------------------
        if DEMParticleGenerator:
            self.x = []
            self.y = []
            self.z = []
            self.m = []
            self.particleRadius = []
            self.compositeParticleIndex = []
            for i in range(SPHGenerator.localNumNodes()):

                xi,yi,zi,mi,Ri = DEMParticleGenerator(SPHGenerator.x[i],
                                                      SPHGenerator.y[i],
                                                      SPHGenerator.z[i],
                                                      SPHGenerator.H[i],
                                                      SPHGenerator.m[i],
                                                      self.particleRadiusFunc(Vector3d(SPHGenerator.x[i],SPHGenerator.y[i],SPHGenerator.z[i])))

                pIDi = [SPHGenerator.globalIDs[i]+1 for entry in Ri]
                self.x.extend(xi)
                self.y.extend(yi)
                self.z.extend(zi)
                self.m.extend(mi)
                self.particleRadius.extend(Ri)
                self.compositeParticleIndex.extend(pIDi)
                
        else:
            self.x = SPHGenerator.x
            self.y = SPHGenerator.y
            self.z = SPHGenerator.z
            self.m = SPHGenerator.m
            self.particleRadius = [self.particleRadiusFunc(Vector3d(SPHGenerator.x[i],SPHGenerator.y[i],SPHGenerator.z[i])) for i in range(len(self.x))]
            self.compositeParticleIndex = SPHGenerator.globalIDs
    
        self.H = [SymTensor3d(1.0/(Rj*hOverR), 0.0, 0.0, 0.0, 1.0/(Rj*hOverR), 0.0, 0.0, 0.0, 1.0/(Rj*hOverR)) for Rj in self.particleRadius]

        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.z, self.m, self.H, self.particleRadius, self.compositeParticleIndex)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localParticleRadius(self, i):
        return self.particleRadius[i]

    def globalParticleRadius(self, i):
        return self._globalValue(i, self.localParticleRadius)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localCompositeParticleIndex(self, i):
        return self.compositeParticleIndex[i]

    def globalCompositeParticleIndex(self, i):
        return self._globalValue(i, self.localCompositeParticleIndex)

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        return Vector3d(self.x[i],self.y[i],self.z[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        return self.H[i]
    
