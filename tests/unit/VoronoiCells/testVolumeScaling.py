#ATS:test(SELF, "--dimension 1", np=1, label="VolumeScaling 1D Cartesian")
#ATS:test(SELF, "--dimension 1 --cartesian False --approxTol 0.1", np=1, label="VolumeScaling 1D Spherical")
#ATS:test(SELF, "--dimension 2", np=1, label="VolumeScaling 2D Cartesian")
#ATS:test(SELF, "--dimension 2 --cartesian False", np=1, label="VolumeScaling 2D RZ")
#ATS:test(SELF, "--dimension 3", np=4, label="VolumeScaling 3D Cartesian")

#-------------------------------------------------------------------------------
# Test that the different volume types give the correct total volume in each
# coordinate system (Cartesian, Spherical, RZ).
#
# Volume is stored as coordinate-plane (unscaled). For a uniform-density
# lattice, the expected total coordinate-plane volume is the domain area/length.
# In non-Cartesian coordinates, mass carries a geometric factor so m/rho gives
# the physical volume; we divide out the geometric factor per node to get the
# expected coordinate-plane volume.
#-------------------------------------------------------------------------------
import os, sys
import numpy as np
from math import *

from Spheral import *
from SpheralTestUtilities import *
import mpi

commandLine(dimension = 1,
            cartesian = True,
            nx = 10,
            x0 = 0.0,
            x1 = 1.5,
            rho0 = 2.0,
            nPerh = 4.01,
            exactTol = 1e-10,
            approxTol = 1e-2,
            skipTol = 1e10)

#-------------------------------------------------------------------------------
# Import the correct version of Spheral
#-------------------------------------------------------------------------------
spherical = False
cylindrical = False
if cartesian:
    exec("from Spheral%id import *" % dimension, globals())
elif dimension == 1:
    from SphericalSpheral import *
    spherical = True
elif dimension == 2:
    from SpheralRZ import *
    cylindrical = True
else:
    raise ValueError("dimension={} for cartesian={} not found".format(dimension, cartesian))

title("VolumeScaling test: dim={}, cartesian={}, spherical={}, cylindrical={}".format(
    dimension, cartesian, spherical, cylindrical))

#-------------------------------------------------------------------------------
# Material and calculation objects
#-------------------------------------------------------------------------------
units = CGuS()
gamma = 5.0/3.0
mu = 1.0
eos = GammaLawGas(gamma, mu, units)

#-------------------------------------------------------------------------------
# Create DataBase and node list
#-------------------------------------------------------------------------------
dataBase = DataBase()
nodes = makeFluidNodeList("nodes", eos,
                          hmin = 1.0e-10,
                          hmax = 1.0e10,
                          nPerh = nPerh,
                          kernelExtent = 2.0)
dataBase.appendNodeList(nodes)
output("nodes")

#-------------------------------------------------------------------------------
# Generate the node distribution
#-------------------------------------------------------------------------------
if dimension == 1:
    if cartesian:
        from DistributeNodes import distributeNodesInRange1d
        distributeNodesInRange1d([(nodes, nx, rho0, (x0, x1))],
                                 nPerh = nPerh)
    else:
        from SortAndDivideRedistributeNodes import distributeNodes1d
        from GenerateSphericalNodeProfile1d import GenerateSphericalNodeProfile1d
        gen = GenerateSphericalNodeProfile1d(nr = nx,
                                             rho = rho0,
                                             rmin = x0,
                                             rmax = x1,
                                             nNodePerh = nPerh)
        distributeNodes1d((nodes, gen))
elif dimension == 2:
    from GenerateNodeDistribution2d import *
    generator = GenerateNodeDistribution2d(distributionType="lattice",
                                           nRadial = nx, nTheta = nx,
                                           xmin = (x0, x0),
                                           xmax = (x1, x1),
                                           rho = rho0,
                                           nNodePerh = nPerh)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes2d
    else:
        from DistributeNodes import distributeNodes2d
    distributeNodes2d((nodes, generator))
else:
    from GenerateNodeDistribution3d import *
    generator = GenerateNodeDistribution3d(distributionType="lattice",
                                           n1 = nx, n2 = nx, n3 = nx,
                                           xmin = (x0, x0, x0),
                                           xmax = (x1, x1, x1),
                                           rho = rho0,
                                           nNodePerh = nPerh)
    if mpi.procs > 1:
        from VoronoiDistributeNodes import distributeNodes3d
    else:
        from DistributeNodes import distributeNodes3d
    distributeNodes3d((nodes, generator))

numInternal = nodes.numInternalNodes
output("numInternal")

#-------------------------------------------------------------------------------
# Compute expected volumes.
# Volume is stored as coordinate-plane (unscaled).
# expectedCoordVol: coordinate-plane domain volume
# expectedPhysVol:  physical volume = totalMass / rho0
# We check that Σ V_i ≈ expectedCoordVol  (coordinate plane)
# and that Σ V_i * geomFactor_i ≈ expectedPhysVol  (physical)
#-------------------------------------------------------------------------------
mass = nodes.mass()
totalMass = mpi.allreduce(sum(mass.internalValues()), mpi.SUM)
expectedPhysVol = totalMass / rho0

if dimension == 1:
    expectedCoordVol = x1 - x0
elif dimension == 2:
    expectedCoordVol = (x1 - x0) * (x1 - x0)
else:
    expectedCoordVol = (x1 - x0)**3

output("totalMass")
output("expectedPhysVol")
output("expectedCoordVol")

#-------------------------------------------------------------------------------
# Boundary conditions
# Cartesian: reflecting on all domain walls
# Spherical: SphericalOriginBoundary + reflecting at outer wall
# RZ: AxisBoundaryRZ + reflecting on non-axis walls
#-------------------------------------------------------------------------------
boundaries = []
if dimension == 1:
    if spherical:
        boundaries.append(SphericalOriginBoundary())
        p1 = Plane(Vector(x1), Vector(-1.0))
        boundaries.append(ReflectingBoundary(p1))
    else:
        p0 = Plane(Vector(x0), Vector( 1.0))
        p1 = Plane(Vector(x1), Vector(-1.0))
        boundaries.append(ReflectingBoundary(p0))
        boundaries.append(ReflectingBoundary(p1))
elif dimension == 2:
    if cylindrical:
        boundaries.append(AxisBoundaryRZ(0.1))
        for d in range(dimension):
            for n, lim in enumerate([x0, x1]):
                if d == 1 and n == 0:
                    continue  # axis BC handles y=0
                point = Vector.zero
                point[d] = lim
                normal = Vector.zero
                normal[d] = 1.0 if n == 0 else -1.0
                boundaries.append(ReflectingBoundary(Plane(point, normal)))
    else:
        for d in range(dimension):
            for n, lim in enumerate([x0, x1]):
                point = Vector.zero
                point[d] = lim
                normal = Vector.zero
                normal[d] = 1.0 if n == 0 else -1.0
                boundaries.append(ReflectingBoundary(Plane(point, normal)))
else:
    for d in range(dimension):
        for n, lim in enumerate([x0, x1]):
            point = Vector.zero
            point[d] = lim
            normal = Vector.zero
            normal[d] = 1.0 if n == 0 else -1.0
            boundaries.append(ReflectingBoundary(Plane(point, normal)))

# If we're parallel we need a distributed boundary
if mpi.procs > 1:
    boundaries.append(TreeDistributedBoundary.instance())

# Get the boundaries ordered coorectly
sortedbcs = sorted(boundaries, key=lambda bc: bc.priority)
boundaries = sortedbcs

output("boundaries")

#-------------------------------------------------------------------------------
# Construct kernels and iterate H to steady state
#-------------------------------------------------------------------------------
if spherical:
    kernel = TableKernel(WendlandC2Kernel(), 1000)
    dimKernel = TableKernel1d(WendlandC2Kernel1d(), 1000)
else:
    kernel = TableKernel(WendlandC2Kernel(), 1000)
    dimKernel = kernel

method = SPHSmoothingScale(IdealH, dimKernel)
iterateIdealH(dataBase,
              [method],
              boundaries,
              100,
              1.e-4)
dataBase.updateConnectivityMap(True, False)

#-------------------------------------------------------------------------------
# Define which volume types to test
#-------------------------------------------------------------------------------
volumeTypes = [(MassOverDensity, exactTol),
               (SumVolume, approxTol),
               (VoronoiVolume, exactTol),
               (HVolume, skipTol),
               (HullVolume, skipTol)]

volumeTypeNames = {MassOverDensity: "MassOverDensity",
                   SumVolume:       "SumVolume",
                   VoronoiVolume:   "VoronoiVolume",
                   HullVolume:      "HullVolume",
                   HVolume:         "HVolume"}

#-------------------------------------------------------------------------------
# For each volume type, create VolumeUpdate (or VoronoiCells), initialize,
# and check the total volume
#-------------------------------------------------------------------------------
checksum = 0

for vt, tol in volumeTypes:
    vtName = volumeTypeNames.get(vt, str(vt))

    # Create the volume package
    if vt == VoronoiVolume:
        volPkg = VoronoiCells(volumeType = VoronoiVolume,
                              W = dimKernel,
                              updateInStep = True,
                              updateInFinalize = False)
    else:
        volPkg = VolumeUpdate(volumeType = vt,
                              W = dimKernel,
                              updateInStep = True,
                              updateInFinalize = False)

    # Attach boundary conditions
    for bc in boundaries:
        volPkg.appendBoundary(bc)

    # Add faceted boundaries for Voronoi (needed to get correct boundary volumes)
    if vt == VoronoiVolume:
        if dimension == 1:
            points = vector_of_Vector([Vector(x0), Vector(x1)])
            facetedBoundary = Box1d(points)
        elif dimension == 2:
            points = vector_of_Vector([Vector(x0, x0), Vector(x0, x1),
                                       Vector(x1, x0), Vector(x1, x1)])
            facetedBoundary = Polygon(points)
        else:
            points = vector_of_Vector([Vector(x0, x0, x0), Vector(x0, x0, x1),
                                       Vector(x0, x1, x0), Vector(x0, x1, x1),
                                       Vector(x1, x0, x0), Vector(x1, x0, x1),
                                       Vector(x1, x1, x0), Vector(x1, x1, x1)])
            facetedBoundary = Polyhedron(points)
        volPkg.addFacetedBoundary(facetedBoundary)

    # Standard Physics initialization
    packages = [volPkg]
    state = State(dataBase, packages)
    derivs = StateDerivatives(dataBase, packages)
    for p in packages:
        p.initializeProblemStartup(dataBase)
        p.registerState(dataBase, state)
        p.registerDerivatives(dataBase, derivs)
    for p in packages:
        p.initializeProblemStartupDependencies(dataBase, state, derivs)
    for p in packages:
        p.preStepInitialize(dataBase, state, derivs)

    # Read out the volume from state and compute both sums
    volume = state.scalarFields(HydroFieldNames.volume)
    pos = nodes.positions()
    localCoordVol = 0.0
    localPhysVol = 0.0
    for k in range(volume.numFields):
        for i in range(volume[k].numInternalElements):
            Vi = volume[k][i]
            localCoordVol += Vi
            if cylindrical:
                ri = abs(pos[i].y)
                localPhysVol += Vi * 2.0 * pi * ri
            elif spherical:
                ri = abs(pos[i].x)
                localPhysVol += Vi * ri * ri
            else:
                localPhysVol += Vi
    totalCoordVol = mpi.allreduce(localCoordVol, mpi.SUM)
    totalPhysVol = mpi.allreduce(localPhysVol, mpi.SUM)

    # Report both comparisons
    coordErr = abs(totalCoordVol - expectedCoordVol) / expectedCoordVol if expectedCoordVol > 0 else abs(totalCoordVol)
    physErr = abs(totalPhysVol - expectedPhysVol) / expectedPhysVol if expectedPhysVol > 0 else abs(totalPhysVol)
    status = "SKIP" if tol >= skipTol else "PASS" if coordErr < tol else "FAIL"
    print("{}: {} coordVol={:.6e} expected={:.6e} relErr={:.4e}  vol3d={:.6e} expected={:.6e} relErr={:.4e}".format(
        status, vtName, totalCoordVol, expectedCoordVol, coordErr, totalPhysVol, expectedPhysVol, physErr))

    if status == "FAIL":
        checksum += 1

        # Print per-node diagnostics for failures
        if numInternal <= 50:
            for i in range(numInternal):
                print("  node {}: pos={} vol={}  m={} rho={}".format(
                    i, pos[i], volume[0][i], mass[i], nodes.massDensity()[i]))

#-------------------------------------------------------------------------------
# Final result
#-------------------------------------------------------------------------------
if checksum > 0:
    raise RuntimeError("{} volume type(s) failed the total volume check".format(checksum))
else:
    print("All volume types passed the total volume check")
