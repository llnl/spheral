#-------------------------------------------------------------------------------
# RejecterBase
#
# This is the base class for Spheral rejecters, which are used to reject nodes
# generators during problem setup.
#
# Rejecter should minimally provide the "accept" method accepting positions as
# tuples of (x,y) [2D] or (x,y,z) [3D], return True/False as to whether the
# given position should be accepted by generator.
#-------------------------------------------------------------------------------
import mpi
from SpheralCompiledPackages import Vector2d, Vector3d

class RejecterBase:

    # In older rejecter would use __call__ to explicitly go
    # through the positions and build new return lists.  We provide that
    # functionality here for backwards compatibility with our newer "accept"
    # style
    def __call__(self, *args):
        assert len(args) in (4, 5), "RejecterBase __call__error: unable to determine dimensionality"
        threeD = len(args) == 5
        n0 = len(args[0])
        for input_arr in args:
            assert len(input_arr) == n0
        
        # Extract the actual arguments
        if threeD:
            x0, y0, z0, m0, H0 = args
        else:
            x0, y0, m0, H0 = args

        # Check if we've been handed serial info (same on all processors), or if
        # the points are already unique per process
        n0test = mpi.allreduce(n0, mpi.MIN) == n0
        n0test = mpi.allreduce(n0test, mpi.MIN)
        if n0test:
            if threeD:
                pos0 = Vector3d(x0[0], y0[0], z0[0])
            else:
                pos0 = Vector2d(x0[0], y0[0])
            pos0test = mpi.allreduce(pos0, mpi.MIN) == pos0
            pos0test = mpi.allreduce(pos0test, mpi.MIN)
            serial = pos0test
        else:
            serial = False

        # If we're serial, we take advantage of any available parallelism to split
        # up the containment testing.  The following algorithm is borrowed 
        # from NodeGeneratorBase to divvy up the ID range.
        if serial:
            ndomain0 = n0//mpi.procs
            remainder = n0 % mpi.procs
            assert remainder < mpi.procs
            ndomain = ndomain0
            if mpi.rank < remainder:
                ndomain += 1
            imin = mpi.rank*ndomain0 + min(mpi.rank, remainder)
            imax = imin + ndomain
        else:
            imin, imax = 0, n0

        # Find the indices of the points we're keeping
        if threeD:
            localIndices = [i for i in range(imin, imax) if self.accept(x0[i], y0[i], z0[i])]
        else:
            localIndices = [i for i in range(imin, imax) if self.accept(x0[i], y0[i])]

        # Cull the values
        x1, y1, z1, m1, H1 = [], [], [], [], []
        for iproc in range(mpi.procs):
            otherIndices = mpi.bcast(localIndices, iproc)
            x1 += [x0[i] for i in otherIndices]
            y1 += [y0[i] for i in otherIndices]
            m1 += [m0[i] for i in otherIndices]
            H1 += [H0[i] for i in otherIndices]
            if threeD:
                z1 += [z0[i] for i in otherIndices]

        # That's it
        if threeD:
            return x1, y1, z1, m1, H1
        else:
            return x1, y1, m1, H1
