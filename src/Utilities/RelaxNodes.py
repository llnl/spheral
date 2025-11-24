import numpy as np
import mpi, time

from SpheralCompiledPackages import *
from spheralDimensions import spheralDimensions
dims = spheralDimensions()

def _RelaxNodesFactory(dimension, cartesian = True):
    assert dimension >= 1 and dimension <= 3

    Base = eval(f"Physics{dimension}d")
    class RelaxNodesImpl(Base):
        """Adjust the node distribution during iterateIdealH"""
        def __init__(self, kernel, nPerh = 4.01,
                     neighborInfo = False,
                     # Options for separating nodes
                     minDist = 0.2, # normalized, so between 0 and 1
                     separateIt = 10, # max num separating steps
                     separateSub = 1, # num sub iterations
                     # Options for relaxing nodes
                     relaxIt = 5,      # max num relaxation steps
                     relaxAlpha = 0.4, # fraction of distance to centroid to move each step
                     # Options for density match
                     densityFunc = None,  # density to match, input is DataBase, output is FieldList of densities
                     densityAlpha = 0.25, # fraction of calculated dist to move particle
                     densityIt = 0,       # num density iterations 
                     densitySub = 4,      # num subcycles for each density it
                     densityTol = 1e-8):  # tolerance for pushing density around
            Base.__init__(self)
            # General args
            self.dimension = dimension
            self.cartesian = cartesian
            self.kernel = kernel
            self.nPerh = nPerh
            self.neighborInfo = neighborInfo

            # Separating
            self.minDist = minDist
            self.separateIt = separateIt
            self.separateSub = separateSub
            assert self.minDist >= 0 and self.minDist <= 1

            # Relaxation
            self.relaxAlpha = relaxAlpha
            self.relaxIt = relaxIt
            assert self.relaxAlpha >= 0 and self.relaxAlpha <= 1

            # Density match
            self.densityFunc = densityFunc
            self.densityAlpha = densityAlpha
            self.densityIt = densityIt
            self.densitySub = densitySub
            self.densityTol = densityTol
            
            # Other stuff
            self.step = 0
            self.VectorType = eval(f"Vector{self.dimension}d")

            # Queue of work
            self.separateEnd = separateIt
            self.relaxEnd = self.separateEnd + relaxIt
            self.densityEnd = self.relaxEnd + densityIt
            
            return

        #---------------------------------------------------------------------------
        # Physics methods
        #---------------------------------------------------------------------------
        # Required methods in Physics that we don't use here
        def evaluateDerivatives(self, t, dt, dataBase, state, derivs):
            return
        def dt(self, dataBase, state, derivs, t):
            return pair_double_string(1e80, self.label())
        def registerState(self, dataBase, state):
            return
        def registerDerivatives(self, dataBase, derivs):
            return
        def label(self):
            return f"RelaxNodes{self.dimension}d"

        # Method that isn't used in iterateIdealH and thus will throw an error
        # if this package is accidentally included elsewhere
        def initializeProblemStartup(self, dataBase):
            raise RuntimeError("RelaxNodes should only be added to iterateIdealH, which (hopefully) doesn't use this function")
        
        # Only methods we actually use
        def preStepInitialize(self, dataBase, state, derivs):
            self.checkSmoothing(dataBase, state, derivs)
            if self.step < self.separateEnd:
                converged = self.separateNodes(dataBase, state, derivs)
                self.step = self.separateEnd if converged else self.step + 1
            elif self.step < self.relaxEnd:
                converged = self.relaxNodes(dataBase, state, derivs, self.relaxAlpha)
                self.step = self.relaxEnd if converged else self.step + 1
            elif self.step < self.densityEnd:
                converged = self.matchDensity(dataBase, state, derivs)
                self.step = self.densityEnd if converged else self.step + 1
            return
        
        def requireVoronoiCells(self):
            return True

        #---------------------------------------------------------------------------
        # Convenience functions
        #---------------------------------------------------------------------------
        def applyBoundaries(self, *fls):
            """Apply boundaries to FieldLists"""
            bcs = self.boundaryConditions
            for bc in bcs:
                for fl in fls:
                    bc.applyFieldListGhostBoundary(fl)
            for bc in bcs:
                bc.finalizeGhostBoundary()
            return

        def getDiameter(self, dataBase, state):
            """Get characteristic diameter of each particle from smoothing length or volume"""
            H = dataBase.fluidHfield
            vol = state.scalarFields(HydroFieldNames.volume)
            nodeLists = dataBase.fluidNodeLists
            self.applyBoundaries(vol, H)
        
            diameter = dataBase.newFluidScalarFieldList()
            dim = self.dimension
            nPerh = self.nPerh
            
            for ni in range(len(nodeLists)):
                for i in range(nodeLists[ni].numNodes):
                    HTr = H[ni][i].Trace()
                    if HTr > 1e-12:
                        diameter[ni][i] = dim / (nPerh * HTr + 1e-12)
                    elif vol[ni][i] > 1e-12:
                        if self.dimension == 1:
                            diameter[ni][i] = vol[ni][i]
                        elif self.dimension == 2:
                            diameter[ni][i] = 2.0 * np.sqrt(vol[ni][i] / np.pi)
                        else:
                            diameter[ni][i] = 2.0 * np.cbrt(0.75 / np.pi * vol[ni][i])
                    else:
                        diameter[ni][i] = 0.0
                        
            return diameter

        def scaleGeom(self, field, dataBase, inverse = False):
            """Remove or add radial scaling to field"""
            if self.cartesian:
                return
            
            if self.dimension == 1:
                raise ValueError("spherical scaling not implemented")
            elif self.dimension == 2:
                pos = dataBase.fluidPosition
                nodeLists = dataBase.fluidNodeLists
                # Need abs for axis boundary
                if inverse:
                    for ni in range(len(nodeLists)):
                        for i in range(nodeLists[ni].numNodes):
                            field[ni][i] /= 2 * np.pi * (abs(pos[ni][i][1]) + 1e-12)
                else:
                    for ni in range(len(nodeLists)):
                        for i in range(nodeLists[ni].numNodes):
                            field[ni][i] *= 2 * np.pi * (abs(pos[ni][i][1]) + 1e-12)
                return
            
            raise ValueError(f"no combination found for dimension={self.dimension} and cartesian={self.cartesian}")
        
        def checkSmoothing(self, dataBase, state, derivs):
            """Check the smoothing length and scale H back to fewer neighbors if needed"""
            if not self.neighborInfo:
                return
            start = time.perf_counter()
            pos = dataBase.fluidPosition
            mass = dataBase.fluidMass
            conn = dataBase.connectivityMap()
            nodeLists = dataBase.fluidNodeLists
            # Hold = state.symTensorFields(HydroFieldNames.H)
            # H = derivs.symTensorFields(SymTensorReplaceBoundedState1d.prefix() + HydroFieldNames.H)
            H = dataBase.fluidHfield
            
            dim = self.dimension
            minDist = self.minDist
            nPerh = self.nPerh

            # Approximate num points should be volume of kernel scaled to nPerh
            if dim == 1:
                goalN = 2.0 * self.nPerh
            elif dim == 2:
                goalN = np.pi * self.nPerh ** 2
            else:
                goalN = 4.0 / 3.0 * np.pi * self.nPerh ** 3
            totN = 0
            minN = 100000000
            maxN = 0
            hmin = 1.0e80
            hmax = 0.0
            htot = 0.0
            for ni in range(len(nodeLists)):
                for i in range(nodeLists[ni].numInternalNodes):
                    currN = conn.numNeighborsForNode(ni, i)
                    totN += currN
                    minN = min(minN, currN)
                    maxN = max(maxN, currN)
                    heig = H[ni][i].eigenValues()
                    h = 1.0 / (abs(H[ni][i].xx) + 1.0e-80)
                    hmin = min(hmin, 1.0 / (abs(heig.maxElement()) + 1.0e-80))
                    hmax = max(hmax, 1.0 / (abs(heig.minElement()) + 1.0e-80))
                    htot += self.dimension / (abs(heig.sumElements()) + 1.0e-80)
            
            numGlob = dataBase.globalNumInternalNodes
            totN = mpi.allreduce(totN, mpi.SUM)
            minN = mpi.allreduce(minN, mpi.MIN)
            maxN = mpi.allreduce(maxN, mpi.MAX)
            if maxN > 3 * goalN:
                meanN = totN / numGlob
                hmin = mpi.allreduce(hmin, mpi.MIN)
                hmax = mpi.allreduce(hmax, mpi.MAX)
                htot = mpi.allreduce(htot, mpi.SUM)
                hmean = htot / numGlob
                print(f"  Neighbor info: (mean, min, max) per node = ({meanN:.6g}, {minN}, {maxN})\n              h: (mean, min, max) = ({hmean:.4e}, {hmin:.4e}, {hmax:.4e})")
            return
        
        #---------------------------------------------------------------------------
        # Separate nodes that are within some fraction of the diameter
        #---------------------------------------------------------------------------
        def getSeparation(self, xij, xijMag, dist):
            """Get the distance vector for a new node separation"""
            if xijMag / dist < 1e-12:
                dirij = self.VectorType.zero
                for d in range(self.dimension):
                    dirij[d] = np.random.uniform(-1, 1)
                dirij = dirij.unitVector()
                return dist * dirij
            sij = xij.unitVector() * (dist - xijMag)
            assert sij.magnitude() < dist + 1e-8
            return sij

        def getSeparationRand(self, xij, xijMag, dist, hasGhost):
            """Get the distance vector for a random separation"""
            # If there are ghosts, need to be deterministic
            if hasGhost:
                return self.getSeparation(xij, xijMag, dist)

            # Pick a random vector
            dirij = self.VectorType.zero
            for d in range(self.dimension):
                dirij[d] = np.random.uniform(-1, 1)
            dirij = dirij.unitVector()

            # If it's really small, just move in the random direction
            if xijMag / dist < 1e-12:
                sij = dist * dirij
                if sij.magnitude() > dist:
                    raise ValueError(f"sij: {sij}")
                assert sij.magnitude() < dist + 1e-8
                return sij

            # Reverse it if it isn't within 90 degrees of xij
            dot = xij.dot(dirij)
            if xij.dot(dirij) < 0:
                dirij = -dirij
                dot = xij.dot(dirij)

            # Compute the distance along that direction needed to separate points
            disc = dot**2 - (xijMag**2 - dist**2)

            # Shouldn't happen, but just in case
            if disc < 0:
                s = dist
                print("warning: quadratic solve wrong")
            else:
                disc2 = np.sqrt(disc)
                s1 = -dot + disc2
                s2 = -dot - disc2
                s = min(x for x in (s1, s2) if x > 0)
            sij = dirij * s
            assert sij.magnitude() < dist + 1e-8
            return sij
        
        def separateNodes(self, dataBase, state, derivs):
            """
            Move the nodes that are too close together apart a bit so their cells
            become unique and can be used for relaxation
            """
            pos = dataBase.fluidPosition
            mass = dataBase.fluidMass
            H = dataBase.fluidHfield
            conn = dataBase.connectivityMap()
            pairs = conn.nodePairList
            nodeLists = dataBase.fluidNodeLists
            numInternal = [nodeLists[ni].numInternalNodes for ni in range(len(nodeLists))]

            self.applyBoundaries(pos, mass)
            self.scaleGeom(mass, dataBase, True)
            
            dim = self.dimension
            minDist = self.minDist
            nPerh = self.nPerh
            diam = self.getDiameter(dataBase, state)
            diamMax = dataBase.newFluidScalarFieldList()

            for it in range(self.separateSub):
                start = time.perf_counter()
                numFixed = 0
                totFrac = 0.0
                maxFrac = 0.0
                for pair in pairs:
                    ni = pair.i_list
                    nj = pair.j_list
                    i = pair.i_node
                    j = pair.j_node

                    # Average the smoothing parameter and get a characteristic length
                    diamij = 0.5 * (diam[ni][i] + diam[nj][j])
                    minDistij = minDist * max(min(diam[ni][i], diam[nj][j]), 0.01 * diamij)
                    diamMax[ni][i] = max(diamMax[ni][i], minDistij)
                    diamMax[nj][j] = max(diamMax[nj][j], minDistij)

                    # If particles are too close, separate them
                    xi = pos[ni][i]
                    xj = pos[nj][j]
                    xij = xi - xj
                    xijMag = xij.magnitude()
                    if xijMag < minDistij:
                        goalDistij = 1.01 * minDistij # to prevent multiple iterations
                        # sij = self.getSeparation(xij, xijMag, goalDistij)
                        hasGhost = [i >= numInternal[ni] or j >= numInternal[nj]]
                        sij = self.getSeparationRand(xij, xijMag, goalDistij, hasGhost)

                        # Preserve the centroid
                        mi = mass[ni][i]
                        mj = mass[nj][j]
                        mijInv = 1.0 / (mi + mj + 1.0e-20)
                        fraci = mj * mijInv
                        fracj = mi * mijInv
                        pos[ni][i] = xi + sij * fraci
                        pos[nj][j] = xj - sij * fracj
                        newDistij = (pos[ni][i] - pos[nj][j]).magnitude()

                        frac = sij.magnitude() / diamij
                        numFixed += 1
                        totFrac += frac
                        maxFrac = max(maxFrac, frac)

                # Make sure particles haven't gotten too close to the axis
                numNeg = 0
                numAxis = 0
                if not self.cartesian:
                    d = 0 if dim == 1 else 1
                    for ni in range(len(nodeLists)):
                        for i in range(nodeLists[ni].numInternalNodes):
                            if pos[ni][i][d] < 0.0:
                                pos[ni][i][d] = abs(pos[ni][i][d])
                                numNeg += 1
                            minDisti = abs(minDist * max(diamMax[ni][i], diam[ni][i]))
                            if pos[ni][i][d] < minDisti:
                                pos[ni][i][d] = minDisti * np.random.uniform(1.01, 1.2)
                                numAxis += 1
                
                # Output stats
                globFixed = mpi.allreduce(numFixed, mpi.SUM)
                globPairs = mpi.allreduce(len(pairs), mpi.SUM)
                globFrac = mpi.allreduce(totFrac, mpi.SUM)
                globMax = mpi.allreduce(maxFrac, mpi.MAX)
                globNeg = mpi.allreduce(numNeg, mpi.SUM)
                globAxis = mpi.allreduce(numAxis, mpi.SUM)
                elapsed = time.perf_counter() - start
                if mpi.rank == 0:
                    message = ""
                    message += f"  Separated {globFixed} / {globPairs} pairs of nodes in {elapsed:.4g} sec"
                    if globFixed > 0:
                        message += f"\n    frac of diam (avg, max, goal) = ({globFrac / (globFixed + 1e-20):.4e}, {globMax:.4e}, {minDist:.4e})"
                    if globNeg > 0 or globAxis > 0:
                        message += f"\n    fixes for axis boundary: num (negative, too close) = ({globNeg}, {globAxis})"
                    print(message)
            # Unscale the mass
            self.scaleGeom(mass, dataBase, False)
            
            # Apply boundaries to affected fields
            self.applyBoundaries(pos)

            return globFixed == 0

        def relaxNodes(self, dataBase, state, derivs, alpha):
            """Relax the nodes using Lloyd's algorithm"""
            assert state.fieldNameRegistered(HydroFieldNames.volume) and state.fieldNameRegistered(HydroFieldNames.cells)
            start = time.perf_counter()
            mass = dataBase.fluidMass
            pos = dataBase.fluidPosition
            H = dataBase.fluidHfield
            rho = dataBase.fluidMassDensity
            vol = state.scalarFields(HydroFieldNames.volume)
            cells = state.facetedVolumeFields(HydroFieldNames.cells)
            nodeLists = dataBase.fluidNodeLists
            flags = state.vector_of_CellFaceFlagFields(HydroFieldNames.cellFaceFlags)

            diam = self.getDiameter(dataBase, state)
            if alpha < 0.0 or alpha > 1.0:
                raise ValueError(f"relax alpha not in bounds: {alpha}")

            sigNum = 0
            totDist = 0.0
            maxDist = 0.0
            for ni in range(len(nodeLists)):
                for i in range(nodeLists[ni].numInternalNodes):
                    # If diam is zero, then there's nothing to relax this time
                    if diam[ni][i] < 1e-12:
                        continue
                    # Check for void on other side and break if found
                    if any(flag.nodeListj == -1 for flag in flags[ni][i]):
                        continue
                    delta = (cells[ni][i].centroid - pos[ni][i]) * alpha
                    pos[ni][i] = pos[ni][i] + delta
                    dist = delta.magnitude() / diam[ni][i]
                    totDist += dist
                    maxDist = max(maxDist, dist)
                    if dist > 1e-8:
                        sigNum += 1

            # Apply boundaries to affected fields
            self.applyBoundaries(pos)

            # Output stats
            elapsed = time.perf_counter() - start
            globSig = mpi.allreduce(sigNum, mpi.SUM)
            globMax = mpi.allreduce(maxDist, mpi.MAX)
            globDist = mpi.allreduce(totDist, mpi.SUM)
            numGlob = dataBase.globalNumInternalNodes
            if mpi.rank == 0:
                message = ""
                message += f"  Relaxed {globSig} / {numGlob} nodes in {elapsed:.4g} sec"
                if globSig > 0:
                    message += f"\n    (mean, max) frac of diam = ({globDist / (numGlob + 1e-20):.4e}, {globMax:.4e})"
                print(message)

            return globSig == 0

        def matchDensity(self, dataBase, state, derivs):
            """Match a given density"""
            if self.densityFunc is None:
                return True
            
            start = time.perf_counter()
            mass = dataBase.fluidMass
            pos = dataBase.fluidPosition
            H = dataBase.fluidHfield
            rho = dataBase.fluidMassDensity
            vol = state.scalarFields(HydroFieldNames.volume)
            flags = state.vector_of_CellFaceFlagFields(HydroFieldNames.cellFaceFlags)
            nodeLists = dataBase.fluidNodeLists
            numNodeLists = dataBase.numFluidNodeLists
            totalNumNodes = dataBase.numInternalNodes

            # Scale the mass if needed
            self.scaleGeom(mass, dataBase, True)

            numMoved = 0.0
            totDist = 0.0
            minDist = 0.0
            maxDist = 0.0
            rhoDiff = dataBase.newFluidScalarFieldList()
            for it in range(self.densitySub):
                # Get density values from the user-supplied function
                rhoV = self.densityFunc(dataBase)
                
                # Get the density difference (what we want to optimize)
                volTot = 0.0
                rhoTot = 0.0
                for ni in range(len(nodeLists)):
                    for i in range(nodeLists[ni].numInternalNodes):
                        if any(flag.nodeListj == -1 for flag in flags[ni][i]):
                            continue
                        rhoDiff[ni][i] = (rho[ni][i] - rhoV[ni][i])
                        volTot += vol[ni][i]
                        rhoTot += rhoV[ni][i] * vol[ni][i]
                volTot = mpi.allreduce(volTot, mpi.SUM)
                rhoTot = mpi.allreduce(rhoTot, mpi.SUM)
                rhoMean = rhoTot / volTot
                
                # Need to apply boundary conditions for the gradient
                self.applyBoundaries(rhoDiff)

                # Get the gradient of the difference
                gradRhoDiff = gradient(rhoDiff, pos, vol, mass, rho, H, self.kernel)

                # Get the diameter for convergence checks and movement limits
                diam = self.getDiameter(dataBase, state)
                alpha = self.densityAlpha
                tol = self.densityTol

                # Move the particles in the negative gradient direction
                for ni in range(len(nodeLists)):
                    for i in range(nodeLists[ni].numInternalNodes):
                        if any(flag.nodeListj == -1 for flag in flags[ni][i]):
                            continue
                        fi = rhoDiff[ni][i]
                        dfi = gradRhoDiff[ni][i]
                        dfiMag = dfi.magnitude()
                        diami = diam[ni][i]

                        # Quit if either is true:
                        # 1. The distance is less than the diameter times tol
                        # 2. The difference is less than the mean density times tol
                        if fi < dfiMag * diami * tol or fi < rhoMean * tol:
                            continue

                        # Don't move more than the diameter
                        disti = alpha * diami * min(1.0, abs(fi / (dfiMag + 1.0e-80)))
                        dfiUnit = dfi.unitVector()

                        pos[ni][i] = pos[ni][i] - dfiUnit * disti
                        numMoved += 1
                        fraci = disti / diami
                        totDist += fraci
                        minDist = min(minDist, fraci)
                        maxDist = max(maxDist, fraci)

                # Update density for the next iteration
                conn = dataBase.connectivityMap()
                densityCalc = eval(f"computeSPHSumMassDensity{self.dimension}d")
                densityCalc(conn, self.kernel, True, pos, mass, H, rho)
            
            # Apply boundaries to affected fields
            self.applyBoundaries(pos, rho)

            # Unscale the mass
            self.scaleGeom(mass, dataBase, False)
            
            # Print info
            globMoved = int(mpi.allreduce(numMoved, mpi.SUM) / self.densitySub)
            globDist = mpi.allreduce(totDist, mpi.SUM) / self.densitySub
            globMin = mpi.allreduce(minDist, mpi.MIN)
            globMax = mpi.allreduce(maxDist, mpi.MAX)
            numGlob = dataBase.globalNumInternalNodes
            elapsed = time.perf_counter() - start
            if mpi.rank == 0:
                print(f"  Density match (sub={self.densitySub}) moved {globMoved} / {numGlob} nodes in {elapsed:.4g} sec\n    (avg, min, max) frac of diam = ({globDist / (globMoved + 1e-20):.4e}, {globMin:.4e}, {globMax:.4e})")
                    
            return globMoved == 0
        
    return RelaxNodesImpl

#-------------------------------------------------------------------------------
# Create the dimension specific particle relaxers that are actually used
#-------------------------------------------------------------------------------
for dim in dims:
    exec(f"RelaxNodes{dim}d = _RelaxNodesFactory({dim})")

if 1 in dims:
    SphericalRelaxNodes = _RelaxNodesFactory(1, False)

if 2 in dims:
    RelaxNodesRZ = _RelaxNodesFactory(2, False)
