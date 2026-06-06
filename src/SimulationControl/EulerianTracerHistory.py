#-------------------------------------------------------------------------------
# This class defines an Eulerian tracer diagnostic for a set of position.
#-------------------------------------------------------------------------------
import os
import mpi
import Spheral
from NodeHistory import NodeHistory

class EulerianTracerHistory(Spheral.RestartableObject):

    # Constructor.
    def __init__(self,
                 geometry,
                 position,
                 samplefunc,
                 W,
                 db,
                 filename,
                 header = None,
                 labels = None,
                 initializefunc = None,
                 weightfunc = None,
                 ):
        self.restart = Spheral.RestartableObject(self)

        # Set up our internal data.
        self.geometry = geometry
        self.position = position
        self.samplefunc = samplefunc
        self.initializefunc = initializefunc
        self.weightfunc = weightfunc
        self.W = W
        self.db = db
        self.filename = filename
        self.labels = labels
        self.cycleHistory = []
        self.timeHistory = []
        self.sampleHistory = []

        # Open the history file.
        self.file = None
        if mpi.rank == 0:
            if os.path.exists(self.filename):
                os.remove(self.filename)
            self.file = open(self.filename, "w")
            assert self.file is not None

            # Write the optional header string.
            if header:
                self.file.write(header + "\n")

            # Write the optional label line
            if labels:
                self.file.write(("# " + ((len(labels) + 2)*'"%20s" ') + "\n") % (("cycle", "time") + labels))

        return

    # This is the method you add to periodic work.
    def sample(self, cycle, ttime, dt):

        # Import the geometry appropriate Spheral types.
        assert self.geometry in ("1d", "2d", "3d", "RZ")
        exec("from Spheral%s import *" % self.geometry)

        # Do we need to initialize anything?
        if self.initializefunc:
            self.initializefunc()

        # How many sample values are we going for?
        for nodeListi, nodeList in enumerate(self.db.fluidNodeLists):
            if nodeList.numNodes > 0:
                nvals = len(self.samplefunc(nodeListi, 0))
        assert nvals > 0

        # Prepare empty slots in the history.
        self.cycleHistory.append(cycle)
        self.timeHistory.append(ttime)
        self.sampleHistory.append([0.0]*nvals)
        Wsum = 0.0

        # Grab position and H FieldLists.
        positions = self.db.globalPosition
        H = self.db.globalHfield
        Hmin = 1e60*SymTensor.one       # Since H is in inverse length, need a big number

        # Prepare the Neighbor information for sampling at this pos, and walk the neighbors.
        numNodeLists = self.db.numFluidNodeLists
        masterLists = vector_of_vector_of_int()
        coarseNeighbors = vector_of_vector_of_int()
        refineNeighbors = vector_of_vector_of_int()
        self.db.setMasterNodeLists(self.position, Hmin, masterLists, coarseNeighbors)
        assert len(coarseNeighbors) == numNodeLists
        self.db.setRefineNodeLists(self.position, Hmin, coarseNeighbors, refineNeighbors)
        for nodeListj in range(numNodeLists):
            for j in refineNeighbors[nodeListj]:

                # Compute the weighting for this position.
                posj = positions(nodeListj, j)
                Hj = H(nodeListj, j)
                Wj = self.W.kernelValue((Hj*(posj - self.position)).magnitude(), 1.0)**2
                if self.weightfunc:
                    Wj *= self.weightfunc(posj)
                Wsum += Wj

                # Use the user supplied method to extract the field values for this (nodeList, index)
                fieldvals = self.samplefunc(nodeListj, j)
                assert len(fieldvals) == nvals

                # Increment the sampled values for this position.
                for i in range(nvals):
                    self.sampleHistory[-1][i] += Wj*fieldvals[i]

        # Normalize the measurements.
        Wsum = max(1.0e-10, mpi.allreduce(Wsum, mpi.SUM))
        for i in range(nvals):
            self.sampleHistory[-1][i] = mpi.allreduce(self.sampleHistory[-1][i], mpi.SUM)/Wsum

        # Update the history file.
        if mpi.rank == 0:
            assert self.file is not None
            samplestr = ""
            for x in self.sampleHistory[-1]:
                samplestr += str(x) + " "
            self.file.write("%i \t %g \t %s\n" % (cycle, ttime, samplestr))
            self.file.flush()

        return

    # Recreate the output file, flushing our full history to it.
    def flushHistory(self):
        if mpi.rank == 0:
            assert self.file is not None
            n = len(self.cycleHistory)
            assert len(self.timeHistory) == n
            assert len(self.sampleHistory) == n
            if mpi.rank == 0:
                for i in range(n):
                    samplestr = ""
                    for x in self.sampleHistory[i]:
                        samplestr += str(x) + " "
                    self.file.write("%i \t %g \t %s\n" % (self.cycleHistory[i],
                                                          self.timeHistory[i],
                                                          samplestr))
            self.file.flush()
        return
                
    # Label for restart.
    def label(self):
        return "EulerianTracerHistory"

    # Write restart.
    def dumpState(self, file, path):
        file.writeObject(self.filename, path + "/filename")
        file.writeObject(self.cycleHistory, path + "/cycleHistory")
        file.writeObject(self.timeHistory, path + "/timeHistory")
        file.writeObject(self.sampleHistory, path + "/sampleHistory")
        return

    # Read restart.
    def restoreState(self, file, path):
        self.filename = file.readObject(path + "/filename")
        self.cycleHistory = file.readObject(path + "/cycleHistory")
        self.timeHistory = file.readObject(path + "/timeHistory")
        self.sampleHistory = file.readObject(path + "/sampleHistory")
        self.flushHistory()
        return
