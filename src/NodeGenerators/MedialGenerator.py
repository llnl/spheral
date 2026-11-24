from math import *
import os
import mpi
import random

from NodeGeneratorBase import *
from SpheralCompiledPackages import *

from sobol import i4_sobol
from centroidalRelaxNodes import centroidalRelaxNodes

#-------------------------------------------------------------------------------
# Base MedialGenerator.  This implements the generic 2D/3D algorithm.
#-------------------------------------------------------------------------------
class MedialGeneratorBase(NodeGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 ndim,
                 n,
                 rho,
                 boundary,
                 gradrho,
                 holes,
                 centroidFrac,
                 maxIterations,
                 fracTol,
                 tessellationBaseDir,
                 tessellationFileName,
                 nNodePerh,
                 randomseed,
                 maxNodesPerDomain,
                 seedPositions,
                 enforceConstantMassPoints,
                 cacheFileName):

        assert ndim in (2,3)
        assert n > 0

        # Load our handy aliases.
        if ndim == 2:
            import Spheral2d as sph
        else:
            import Spheral3d as sph

        # Did we get passed a function or a constant for the density?
        if type(rho) in (float, int):
            def rhofunc(posi):
                return rho
            rhomax = rho
        else:
            rhofunc = rho
            rhomax = None
        self.rhofunc = rhofunc

        # Some useful geometry.
        box = boundary.xmax - boundary.xmin
        length = box.maxElement()
        boundvol = boundary.volume
        for hole in holes:
            boundvol -= hole.volume
        if boundvol <= 0.0:
            # The holes were not entirely contained in the bounding volume, so we punt.
            boundvol = 0.5*boundary.volume
        boxvol = 1.0
        for idim in range(ndim):
            boxvol *= box[idim]
        fracOccupied = min(1.0, boxvol/boundvol)
        assert fracOccupied > 0.0 and fracOccupied <= 1.0

        # If there is an pre-existing cache file, load it instead of doing all the work.
        if not self.restoreState(cacheFileName):

            # Create a temporary NodeList we'll use store and update positions.
            eos = sph.GammaLawGasMKS(2.0, 2.0)
            WT = sph.TableKernel(sph.NBSplineKernel(7), 1000)
            hmax = 2.0*(boundvol/pi*n)**(1.0/ndim)
            nodes = sph.makeFluidNodeList("tmp generator nodes", 
                                          eos,
                                          hmin = 1e-10,
                                          hmax = hmax,
                                          kernelExtent = WT.kernelExtent,
                                          hminratio = 1.0,
                                          nPerh = nNodePerh,
                                          topGridCellSize = 2.0*WT.kernelExtent*hmax)
        
            # Make a first pass looking for the maximum density (roughly).
            pos = nodes.positions()
            mass = nodes.mass()
            rhof = nodes.massDensity()
            H = nodes.Hfield()
            imin, imax = self.globalIDRange(n)
            nlocal = imax - imin
            nodes.numInternalNodes = nlocal
        
            # If the user provided the starting or seed positions, use 'em.
            if seedPositions is not None:
                hi = min(hmax, 2.0 * (boundvol/(pi*n))**(1.0/ndim))
                assert hi > 0.0
                nlocal = len(seedPositions)
                assert mpi.allreduce(nlocal, mpi.SUM) == n
                nodes.numInternalNodes = nlocal
                for i in range(nlocal):
                    pos[i] = seedPositions[i]
                    rhoi = rhofunc(pos[i])
                    rhof[i] = rhoi
                    mass[i] = rhoi * boundvol/n  # Not actually correct, but mass will be updated in centroidalRelaxNodes
                    H[i] = sph.SymTensor.one / hi
        
            else:
                # If necessary probe for a maximum density statistically.
                rangen = random.Random(randomseed + mpi.rank)
                if not rhomax:
                    rhomax = 0.0
                    nglobal = 0
                    while nglobal < n:
                        p = boundary.xmin + length*sph.Vector(rangen.random(), rangen.random(), rangen.random())
                        use = boundary.contains(p, False)
                        if use:
                            ihole = 0
                            while use and ihole < len(holes):
                                use = not holes[ihole].contains(p, True)
                                ihole += 1
                        if use:
                            rhomax = max(rhomax, rhofunc(p))
                            i = 1
                        else:
                            i = 0
                        nglobal += mpi.allreduce(i, mpi.SUM)
                    rhomax = mpi.allreduce(rhomax, mpi.MAX)
                print("MedialGenerator: selected a maximum density of ", rhomax)
            
                # It's a bit tricky to properly use the Sobol sequence in parallel.  We handle this by searching for the lowest
                # seeds that give us the desired number of points.
                seeds = []
                seed = 0
                while mpi.allreduce(len(seeds), mpi.SUM) < n:
                    localseed = seed + mpi.rank
                    [coords, newseed] = i4_sobol(ndim, localseed)
                    p = boundary.xmin + length*sph.Vector(*tuple(coords))
                    use = boundary.contains(p, False)
                    if use:
                        ihole = 0
                        while use and ihole < len(holes):
                            use = not holes[ihole].contains(p, True)
                            ihole += 1
                    if use:
                        rhoi = rhofunc(p)
                        if rangen.random() < rhoi/rhomax:
                            seeds.append(localseed)
                    seed += mpi.procs
            
                # Drop the highest value seeds to ensure we have the correct number of total points.
                nglobal = mpi.allreduce(len(seeds), mpi.SUM)
                assert n + mpi.procs >= nglobal
                seeds.sort()
                seeds = [-1] + seeds
                while mpi.allreduce(len(seeds), mpi.SUM) > n + mpi.procs:
                    maxseed = mpi.allreduce(seeds[-1], mpi.MAX)
                    assert maxseed > -1
                    if seeds[-1] == maxseed:
                        seeds = seeds[:-1]
                seeds = seeds[1:]
            
                # Load balance the number of seeds per domain.
                if len(seeds) > nlocal:
                    extraseeds = seeds[nlocal:]
                else:
                    extraseeds = []
                extraseeds = mpi.allreduce(extraseeds, mpi.SUM)
                seeds = seeds[:nlocal]
                for iproc in range(mpi.procs):
                    ngrab = max(0, nlocal - len(seeds))
                    ntaken = mpi.bcast(ngrab, root=iproc)
                    if mpi.rank == iproc:
                        seeds += extraseeds[:ngrab]
                    extraseeds = extraseeds[ntaken:]
                assert len(extraseeds) == 0
                assert len(seeds) == nlocal
                assert mpi.allreduce(len(seeds), mpi.SUM) == n
            
                # Initialize the desired number of generators in the boundary using the Sobol sequence.
                hi = min(hmax, 2.0 * (boundvol/(pi*n))**(1.0/ndim))
                assert hi > 0.0
                for i, seed in enumerate(seeds):
                    [coords, newseed] = i4_sobol(ndim, seed)
                    p = boundary.xmin + length*sph.Vector(*tuple(coords))
                    rhoi = rhofunc(p)
                    pos[i] = p
                    rhof[i] = rhoi
                    mass[i] = rhoi * boundvol/n  # Not actually correct, but mass will be updated in centroidalRelaxNodes
                    H[i] = sph.SymTensor.one / hi
        
                # Each domain has independently generated the correct number of points, but they are randomly distributed.
                # Before going further it's useful to try and spatially collect the points by domain.
                # We'll use the Spheral Peano-Hilbert space filling curve implementation to do this.
                if mpi.procs > 1:
                    db = sph.DataBase()
                    db.appendNodeList(nodes)
                    maxNodes = max(maxNodesPerDomain, 2*n/mpi.procs)
                    redistributor = sph.PeanoHilbertOrderRedistributeNodes()
                    redistributor.redistributeNodes(db)
        
            # If we're in parallel we need the parallel boundary.
            if mpi.procs > 1:
                boundaries = [sph.TreeDistributedBoundary.instance()]
            else:
                boundaries = []
        
            # Iterate the points toward centroidal relaxation.
            vol, surfacePoint = centroidalRelaxNodes([(nodes, boundary, holes)],
                                                     W = WT,
                                                     rho = rhofunc,
                                                     gradrho = gradrho,
                                                     boundaries = boundaries,
                                                     fracTol = fracTol,
                                                     centroidFrac = centroidFrac,
                                                     maxIterations = maxIterations,
                                                     tessellationBaseDir = tessellationBaseDir,
                                                     tessellationFileName = tessellationFileName)
        
            # Store the values the descendent generators will need.
            self.vol, self.surface, self.pos, self.m, self.H = [], [], [], [], []
            for i in range(nodes.numInternalNodes):
                self.vol.append(vol(0,i))
                self.surface.append(surfacePoint(0,i))
                self.pos.append(sph.Vector(pos[i]))
                self.m.append(vol(0,i) * rhofunc(pos[i]))
                self.H.append(sph.SymTensor(H[i]))
            assert mpi.allreduce(len(self.vol), mpi.SUM) == n
            assert mpi.allreduce(len(self.surface), mpi.SUM) == n
            assert mpi.allreduce(len(self.pos), mpi.SUM) == n
            assert mpi.allreduce(len(self.m), mpi.SUM) == n
            assert mpi.allreduce(len(self.H), mpi.SUM) == n
        
            # If requested, enforce constant mass points.
            if enforceConstantMassPoints:
                msum = mpi.allreduce(sum([0.0] + self.m), mpi.SUM)
                self.m = [msum/n]*len(self.pos)

            # If requested, we can store the state of the generator such that it can be
            # later restored without going through all that work.
            if cacheFileName:
                self.dumpState(cacheFileName)

        return

    #---------------------------------------------------------------------------
    # Try to read the state from a file.
    #---------------------------------------------------------------------------
    def restoreState(self, cacheFileName):
        def readNodeData(f, iproc):
            pos = f.readObject("proc%06i/pos" % iproc)
            m = f.readObject("proc%06i/m" % iproc)
            H = f.readObject("proc%06i/H" % iproc)
            vol = f.readObject("proc%06i/vol" % iproc)
            surface = f.readObject("proc%06i/surface" % iproc)
            return pos, m, H, vol, surface

        if cacheFileName is None:
            return False

        if os.path.splitext(cacheFileName) != ".silo":
            cacheFileName += ".silo"
        result = False
        if mpi.rank == 0:
            result = (cacheFileName and os.path.exists(cacheFileName))
        result = mpi.bcast(result, root=0)
        if result:
            print("Restoring MedialGenerator state from %s" % cacheFileName)
            if mpi.rank == 0:
                f = SiloFileIO(cacheFileName, Read)
                numGeneratingProcs = f.readObject("numGeneratingProcs")
                
                # Decide how to divide the generating domains between our current processes.
                n0 = numGeneratingProcs/mpi.procs
                remainder = numGeneratingProcs % mpi.procs
                for iproc in range(mpi.procs):
                    if iproc >= numGeneratingProcs:
                        imin, imax = 0, 0
                    else:
                        imin = iproc*n0 + min(iproc, remainder)
                        imax = imin + n0
                        if iproc < remainder:
                            imax += 1
                    pos, m, H, vol, surface = [], [], [], [], []
                    for igenproc in range(imin, imax):
                        posi, mi, Hi, voli, surfacei = readNodeData(f, igenproc)
                        pos += posi
                        m += mi
                        H += Hi
                        vol += voli
                        surface += surfacei
                    if iproc == 0:
                        self.pos, self.m, self.H, self.vol, self.surface = pos, m, H, vol, surface
                    else:
                        mpi.send((pos, m, H, vol, surface), dest=iproc)
                f.close()
            else:
                self.pos, self.m, self.H, self.vol, self.surface = mpi.recv(source=0)[0]
        return result

    #---------------------------------------------------------------------------
    # Write the state to a file.
    #---------------------------------------------------------------------------
    def dumpState(self, cacheFileName):
        def writeNodeData(f, iproc, pos, m, H, vol, surface):
            f.writeObject(pos, "proc%06i/pos" % iproc)
            f.writeObject(m, "proc%06i/m" % iproc)
            f.writeObject(H, "proc%06i/H" % iproc)
            f.writeObject(vol, "proc%06i/vol" % iproc)
            f.writeObject(surface, "proc%06i/surface" % iproc)
            return

        if os.path.splitext(cacheFileName) != ".silo":
            cacheFileName += ".silo"
        if mpi.rank == 0:
            dire = os.path.dirname(cacheFileName)
            if dire and not os.path.exists(dire):
                os.makedirs(dire)
            f = SiloFileIO(cacheFileName, Create)
            f.writeObject(mpi.procs, "numGeneratingProcs")
            writeNodeData(f, 0, self.pos, self.m, self.H, self.vol, self.surface)
            for iproc in range(1, mpi.procs):
                pos, m, H, vol, surface = mpi.recv(source=iproc)[0]
                writeNodeData(f, iproc, pos, m, H, vol, surface)
            f.close()
        else:
            mpi.send((self.pos, self.m, self.H, self.vol, self.surface), dest=0)
        mpi.barrier()
        return

#-------------------------------------------------------------------------------
# 2D Generator.  Seeds positions within the given boundary using the Sobol 
# sequence to create the seeds.
# Note as currently implemented this code does not necessarily generate 
# identical output for different numbers of processors!
#-------------------------------------------------------------------------------
class MedialGenerator2d(MedialGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 n,
                 rho,
                 boundary,
                 gradrho = None,
                 holes = [],
                 centroidFrac = 1.0,
                 maxIterations = 100,
                 fracTol = 1.0e-3,
                 tessellationBaseDir = ".",
                 tessellationFileName = None,
                 nNodePerh = 2.01,
                 offset = (0.0, 0.0),
                 rejecter = None,
                 randomseed = 492739149274,
                 maxNodesPerDomain = 1000,
                 seedPositions = None,
                 enforceConstantMassPoints = False,
                 cacheFileName = None):

        # The base generator does most of the work.
        MedialGeneratorBase.__init__(self,
                                     ndim = 2,
                                     n = n,
                                     rho = rho,
                                     boundary = boundary,
                                     gradrho = gradrho,
                                     holes = holes,
                                     centroidFrac = centroidFrac,
                                     maxIterations = maxIterations,
                                     fracTol = fracTol,
                                     tessellationBaseDir = tessellationBaseDir,
                                     tessellationFileName = tessellationFileName,
                                     nNodePerh = nNodePerh,
                                     randomseed = randomseed,
                                     maxNodesPerDomain = maxNodesPerDomain,
                                     seedPositions = seedPositions,
                                     enforceConstantMassPoints = enforceConstantMassPoints,
                                     cacheFileName = cacheFileName)

        # Convert to our now regrettable standard coordinate storage for generators.
        self.x = [x.x + offset[0] for x in self.pos]
        self.y = [x.y + offset[1] for x in self.pos]
        del self.pos

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.m, self.H, self.vol, self.surface = rejecter(self.x,
                                                                              self.y,
                                                                              self.m,
                                                                              self.H,
                                                                              self.vol,
                                                                              self.surface)

        # Initialize the base class, taking into account the fact we've already broken
        # up the nodes between domains.
        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.m, self.H, self.vol, self.surface)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector2d(self.x[i], self.y[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.H)
        return self.rhofunc(Vector2d(self.x[i], self.y[i]))
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

#-------------------------------------------------------------------------------
# 3D Generator.  Seeds positions within the given boundary using the Sobol 
# sequence to create the seeds.
# Note as currently implemented this code does not necessarily generate 
# identical output for different numbers of processors!
#-------------------------------------------------------------------------------
class MedialGenerator3d(MedialGeneratorBase):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 n,
                 rho,
                 boundary,
                 gradrho = None,
                 holes = [],
                 maxIterations = 100,
                 centroidFrac = 1.0,
                 fracTol = 1.0e-3,
                 tessellationBaseDir = ".",
                 tessellationFileName = None,
                 nNodePerh = 2.01,
                 offset = (0.0, 0.0, 0.0),
                 rejecter = None,
                 randomseed = 492739149274,
                 maxNodesPerDomain = 1000,
                 seedPositions = None,
                 enforceConstantMassPoints = False,
                 cacheFileName = None):

        # The base generator does most of the work.
        MedialGeneratorBase.__init__(self,
                                     ndim = 3,
                                     n = n,
                                     rho = rho,
                                     boundary = boundary,
                                     gradrho = gradrho,
                                     holes = holes,
                                     centroidFrac = centroidFrac,
                                     maxIterations = maxIterations,
                                     fracTol = fracTol,
                                     tessellationBaseDir = tessellationBaseDir,
                                     tessellationFileName = tessellationFileName,
                                     nNodePerh = nNodePerh,
                                     randomseed = randomseed,
                                     maxNodesPerDomain = maxNodesPerDomain,
                                     seedPositions = seedPositions,
                                     enforceConstantMassPoints = enforceConstantMassPoints,
                                     cacheFileName = cacheFileName)

        # Convert to our now regrettable standard coordinate storage for generators.
        self.x = [x.x + offset[0] for x in self.pos]
        self.y = [x.y + offset[1] for x in self.pos]
        self.z = [x.z + offset[2] for x in self.pos]
        del self.pos

        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.z, self.m, self.H, self.vol, self.surface = rejecter(self.x,
                                                                                      self.y,
                                                                                      self.z,
                                                                                      self.m,
                                                                                      self.H,
                                                                                      self.vol,
                                                                                      self.surface)

        # Initialize the base class, taking into account the fact we've already broken
        # up the nodes between domains.
        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.z, self.m, self.H, self.vol, self.surface)
        return

    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector3d(self.x[i], self.y[i], self.z[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.H)
        return self.rhofunc(Vector3d(self.x[i], self.y[i], self.z[i]))
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]

#-------------------------------------------------------------------------------
# 3D Generator.  Seeds positions using the RPRPS algorithm.
#-------------------------------------------------------------------------------
class MedialSphereGenerator3d(MedialGeneratorBase):
    
    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 n,
                 rho,
                 rmin,
                 rmax,
                 gradrho = None,
                 maxIterations = 100,
                 centroidFrac = 1.0,
                 fracTol = 1.0e-3,
                 tessellationBaseDir = ".",
                 tessellationFileName = None,
                 nNodePerh = 2.01,
                 offset = (0.0, 0.0, 0.0),
                 rejecter = None,
                 randomseed = 492739149274,
                 maxNodesPerDomain = 1000,
                 seedPositions = None,
                 enforceConstantMassPoints = False,
                 cacheFileName = None):
        
        from Spheral3d import vector_of_Vector, Vector
        from GenerateNodeDistribution3d import GenerateIcosahedronMatchingProfile3d
        SphericallyConformalMap = GenerateIcosahedronMatchingProfile3d(n=n,
                                                                       densityProfileMethod=rho,
                                                                       rmin=rmin,
                                                                       rmax=rmax,
                                                                       nNodePerh=nNodePerh,
                                                                       offset=offset,
                                                                       rejecter=rejecter)
        n = len(SphericallyConformalMap.positions)
        # this is a serialized list, so it needs to be broken up to each processor
        # no degeneracies allowed!
        seedPositions = vector_of_Vector()
        nprocs = mpi.procs
        nlocal = int(n/nprocs)
        nrem   = n - nlocal*nprocs
        myn    = nlocal
        if mpi.rank == nprocs-1:
            myn += nrem
        for i in range(myn):
            j = i + mpi.rank*nlocal
            seedPositions.append(Vector(SphericallyConformalMap.positions[j][0],
                                        SphericallyConformalMap.positions[j][1],
                                        SphericallyConformalMap.positions[j][2]))
        
        
        # Now construct boundaries based on rmin and rmax
        bpoints = vector_of_Vector()
        nshell = 300 # fix this later
        for i in range(1,nshell+1):
            h = -1.0+(2.0*(i-1.0)/(nshell-1.0))
            t = acos(h)
            if (i>1 and i<nshell):
                p = (p + 3.8/sqrt(nshell)*1.0/sqrt(1.0-h*h)) % (2.0*pi)
            else:
                p = 0
            bpoints.append(Vector(rmax*sin(t)*cos(p),
                                  rmax*sin(t)*sin(p),
                                  rmax*cos(t)))
        bound1 = Polyhedron(bpoints)

        if (rmin > 0.0):
            for i in range(nshell):
                bpoints[i] *= rmin/rmax
            bound2 = Polyhedron(bpoints)
        
        # The base generator does most of the work.
        MedialGeneratorBase.__init__(self,
                                     ndim = 3,
                                     n = n,
                                     rho = rho,
                                     boundary = bound1,
                                     gradrho = gradrho,
                                     holes = bound2 if (rmin>0.0) else [],
                                     centroidFrac = centroidFrac,
                                     maxIterations = maxIterations,
                                     fracTol = fracTol,
                                     tessellationBaseDir = tessellationBaseDir,
                                     tessellationFileName = tessellationFileName,
                                     nNodePerh = nNodePerh,
                                     randomseed = randomseed,
                                     maxNodesPerDomain = maxNodesPerDomain,
                                     seedPositions = seedPositions,
                                     enforceConstantMassPoints = enforceConstantMassPoints,
                                     cacheFileName = cacheFileName)
    
        # Convert to our now regrettable standard coordinate storage for generators.
        self.x = [x.x + offset[0] for x in self.pos]
        self.y = [x.y + offset[1] for x in self.pos]
        self.z = [x.z + offset[2] for x in self.pos]
        del self.pos
        
        # If the user provided a "rejecter", give it a pass
        # at the nodes.
        if rejecter:
            self.x, self.y, self.z, self.m, self.H, self.vol, self.surface = rejecter(self.x,
                                                                                      self.y,
                                                                                      self.z,
                                                                                      self.m,
                                                                                      self.H,
                                                                                      self.vol,
                                                                                      self.surface)

        # Initialize the base class, taking into account the fact we've already broken
        # up the nodes between domains.
        NodeGeneratorBase.__init__(self, False,
                                   self.x, self.y, self.z, self.m, self.H, self.vol, self.surface)
        return
    
    #---------------------------------------------------------------------------
    # Get the position for the given node index.
    #---------------------------------------------------------------------------
    def localPosition(self, i):
        assert i >= 0 and i < len(self.x)
        assert len(self.x) == len(self.y)
        return Vector3d(self.x[i], self.y[i], self.z[i])
    
    #---------------------------------------------------------------------------
    # Get the mass for the given node index.
    #---------------------------------------------------------------------------
    def localMass(self, i):
        assert i >= 0 and i < len(self.m)
        return self.m[i]
    
    #---------------------------------------------------------------------------
    # Get the mass density for the given node index.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        assert i >= 0 and i < len(self.H)
        return self.rhofunc(Vector3d(self.x[i], self.y[i], self.z[i]))
    
    #---------------------------------------------------------------------------
    # Get the H tensor for the given node index.
    #---------------------------------------------------------------------------
    def localHtensor(self, i):
        assert i >= 0 and i < len(self.H)
        return self.H[i]
