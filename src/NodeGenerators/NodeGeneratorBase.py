from math import *

import mpi
rank = mpi.rank
procs = mpi.procs

from Spheral import Vector1d, SymTensor1d, Vector2d, SymTensor2d, Vector3d, SymTensor3d, rotationMatrix2d, rotationMatrix3d

#-------------------------------------------------------------------------------
# NodeGeneratorBase.
#
# Base class for Spheral++ NodeGenerators.  Implements the default methods
# for getting node info by global IDs and such.
#-------------------------------------------------------------------------------
class NodeGeneratorBase:

    def __init__(self,
                 serialInitialization,
                 *vars):

        if serialInitialization:
            ntot = len(vars[0])
            minGlobalID, maxGlobalID = self.globalIDRange(ntot)
            self.globalIDs = list(range(minGlobalID, maxGlobalID))
            self._cullVars(minGlobalID, maxGlobalID, *vars)

        else:
            # One allgather of the local counts (rather than a bcast per rank)
            # gives each rank the offset of its global IDs.
            counts = mpi.allgather(len(vars[0]))
            offset = sum(counts[:mpi.rank])
            self.globalIDs = list(range(offset, offset + len(vars[0])))

        return

    #---------------------------------------------------------------------------
    # Simple minded method to compute and assign unique global ID ranges to
    # each processor.
    #---------------------------------------------------------------------------
    def globalIDRange(self, ntot):        
        ndomain0 = ntot // procs
        remainder = ntot % procs
        assert remainder < procs
        ndomain = ndomain0
        if rank < remainder:
            ndomain += 1
        minGlobalID = rank*ndomain0 + min(rank, remainder)
        maxGlobalID = minGlobalID + ndomain
        
        assert mpi.allreduce(minGlobalID, mpi.MIN) == 0
        assert mpi.allreduce(maxGlobalID, mpi.MAX) == ntot

        return minGlobalID, maxGlobalID

    #---------------------------------------------------------------------------
    # Cull the given set of variables to the given indicies.
    #---------------------------------------------------------------------------
    def _cullVars(self, minID, maxID, *args):
        for arg in args:
            del arg[maxID:]
            del arg[:minID]
        return

    #---------------------------------------------------------------------------
    # Find the domain that possesses the given globalID.
    #---------------------------------------------------------------------------
    def domainForID(self, i):
        result = -1
        if i in self.globalIDs:
            result = rank
        mpi.allreduce(result, mpi.MAX)
        return result

    #---------------------------------------------------------------------------
    # Get the position for the given global node index.
    #---------------------------------------------------------------------------
    def globalPosition(self, i):
        return self._globalValue(i, self.localPosition)

    #-------------------------------------------------------------------------------
    # Get the mass for the given global node index.
    #-------------------------------------------------------------------------------
    def globalMass(self, i):
        return self._globalValue(i, self.localMass)

    #-------------------------------------------------------------------------------
    # Get the mass density for the given global node index.
    #-------------------------------------------------------------------------------
    def globalMassDensity(self, i):
        return self._globalValue(i, self.localMassDensity)

    #-------------------------------------------------------------------------------
    # Get the H tensor for the given global node index.
    #-------------------------------------------------------------------------------
    def globalHtensor(self, i):
        return self._globalValue(i, self.localHtensor)

    #-------------------------------------------------------------------------------
    # Get the local number of nodes.
    #-------------------------------------------------------------------------------
    def localNumNodes(self):
        return len(self.m)

    #-------------------------------------------------------------------------------
    # Get the total number of nodes.
    #-------------------------------------------------------------------------------
    def globalNumNodes(self):
        return mpi.allreduce(len(self.m), mpi.SUM)

    #---------------------------------------------------------------------------
    # Force the H tensors to be round.
    #---------------------------------------------------------------------------
    def makeHround(self):
        n = len(self.H)
        if n > 0:
            SymTensor = {type(SymTensor2d.zero) : SymTensor2d,
                         type(SymTensor3d.zero) : SymTensor3d}[type(self.H[0])]
            for i in range(n):
                h0 = self.H[i].Inverse().Trace() / SymTensor.nDimensions
                assert h0 > 0.0
                self.H[i] = 1.0/h0 * SymTensor.one
        return

    #---------------------------------------------------------------------------
    # Generic method to get the global value.
    #---------------------------------------------------------------------------
    def _globalValue(self, i, method):
        dom = self.domainForID(i)
        result = None
        if rank == dom:
            j = self.globalIDs[i]
            result = method(j)
        mpi.bcast(result, dom)
        return result

    #---------------------------------------------------------------------------
    # Copy a list to a std::vector.
    #---------------------------------------------------------------------------
    def vectorFromList(self, llist, constructor):
        result = constructor()
        for x in llist:
            result.append(x)
        return result

    #---------------------------------------------------------------------------
    # Default no-op field methods, which can be over-ridden in descendents.
    #---------------------------------------------------------------------------
    def localMassDensity(self, i):
        return 0.0

    def localVelocity(self, i):
        return {type(SymTensor1d.zero) : Vector1d.zero,
                type(SymTensor2d.zero) : Vector2d.zero,
                type(SymTensor3d.zero) : Vector3d.zero}[type(self.H[0])]

#-------------------------------------------------------------------------------
# Helper class for providing default constant density behaviour.
#-------------------------------------------------------------------------------
class ConstantRho:
    def __init__(self, rho0):
        assert type(rho0) is float
        self.rho0 = rho0
        return
    def rho(self, r):
        return self.rho0
    def __call__(self, r):
        return self.rho(r)

#-------------------------------------------------------------------------------
# Helper class for providing default list of densities behaviour.
#-------------------------------------------------------------------------------
class ListRho:
    def __init__(self, rho0):
        assert type(rho0) is list
        self.rho0 = rho0
        return
    def __call__(self, i):
        return self.rho[i]

#-------------------------------------------------------------------------------
# A simple minded node refinement algorithm for a NodeGenerator in 2-D.
#-------------------------------------------------------------------------------
def refineNodes2d(gen,
                  deta = 0.25):
    n = gen.localNumNodes()
    x = []
    y = []
    m = []
    rho = []
    vx = []
    vy = []
    eps = []
    H = []
    extras = {}
    for name in gen.extraFields:
        extras[name] = []

    etas = [Vector2d(-1, -1) * deta,
            Vector2d(-1,  1) * deta,
            Vector2d( 1, -1) * deta,
            Vector2d( 1,  1) * deta]

    for i in range(n):
        ri = gen.localPosition(i)
        mi = gen.localMass(i)
        Hi = gen.localHtensor(i)
        Hinv = Hi.Inverse()
        eigen = Hinv.eigenVectors()
        R = rotationMatrix2d(eigen.eigenVectors.getColumn(0))
        T = SymTensor2d(eigen.eigenValues.x, 0.0,
                        0.0,                 eigen.eigenValues.y)
        T.rotationalTransform(eigen.eigenVectors)
        T = T*R.Inverse()

        mj = 0.25*mi
        Hj = Hi*2.0
        for eta in etas:
            rj = ri + T*eta
            x.append(rj.x)
            y.append(rj.y)
            m.append(mj)
            rho.append(gen.rho[i])
            vx.append(gen.vx[i])
            vy.append(gen.vy[i])
            eps.append(gen.eps[i])
            H.append(Hj)
            for name in gen.extraFields:
                extras[name].append(gen.__dict__[name][i])

    gen.x = x
    gen.y = y
    gen.m = m
    gen.rho = rho
    gen.vx = vx
    gen.vy = vy
    gen.eps = eps
    gen.H = H
    for name in gen.extraFields:
        gen.__dict__[name] = extras[name]

    for f in ([gen.x, gen.y, gen.m, gen.vx, gen.vy, gen.eps, gen.H] +
              [gen.__dict__[x] for x in gen.extraFields]):
        assert len(f) == len(etas)*n
    return

#-------------------------------------------------------------------------------
# A simple minded node refinement algorithm for a NodeGenerator in 3-D.
#-------------------------------------------------------------------------------
def refineNodes3d(gen,
                  deta = 0.25):
    n = gen.localNumNodes()
    x = []
    y = []
    z = []
    m = []
    rho = []
    vx = []
    vy = []
    vz = []
    eps = []
    H = []
    extras = {}
    for name in gen.extraFields:
        extras[name] = []

    etas = [Vector3d(-1, -1, -1) * deta,
            Vector3d(-1, -1,  1) * deta,
            Vector3d(-1,  1, -1) * deta,
            Vector3d(-1,  1,  1) * deta,
            Vector3d( 1, -1, -1) * deta,
            Vector3d( 1, -1,  1) * deta,
            Vector3d( 1,  1, -1) * deta,
            Vector3d( 1,  1,  1) * deta]

    for i in range(n):
        ri = gen.localPosition(i)
        mi = gen.localMass(i)
        Hi = gen.localHtensor(i)
        Hinv = Hi.Inverse()
        eigen = Hinv.eigenVectors()
        R = rotationMatrix3d(eigen.eigenVectors.getColumn(0))
        T = SymTensor2d(eigen.eigenValues.x, 0.0,                 0.0,
                        0.0,                 eigen.eigenValues.y, 0.0,
                        0.0,                 0.0,                 eigen.eigenValues.z)
        T.rotationalTransform(eigen.eigenVectors)
        T = T*R.Inverse()

        mj = mi/8.0
        Hj = Hi*2.0
        for eta in etas:
            rj = ri + T*eta
            x.append(rj.x)
            y.append(rj.y)
            z.append(rj.z)
            m.append(mj)
            rho.append(gen.rho[i])
            vx.append(gen.vx[i])
            vy.append(gen.vy[i])
            vz.append(gen.vz[i])
            eps.append(gen.eps[i])
            H.append(Hj)
            for name in gen.extraFields:
                extras[name].append(gen.__dict__[name][i])

    gen.x = x
    gen.y = y
    gen.z = z
    gen.m = m
    gen.rho = rho
    gen.vx = vx
    gen.vy = vy
    gen.vz = vz
    gen.eps = eps
    gen.H = H
    for name in gen.extraFields:
        gen.__dict__[name] = extras[name]

    for f in ([gen.x, gen.y, gen.z, gen.m, gen.vx, gen.vy, gen.vz, gen.eps, gen.H] +
              [gen.__dict__[x] for x in gen.extraFields]):
        assert len(f) == len(etas)*n
    return
