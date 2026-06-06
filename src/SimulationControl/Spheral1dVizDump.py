#-------------------------------------------------------------------------------
# A 1D analogue to our 2D and 3D Voronoi silo dumps.  This just outputs columns
# of numbers in a column oriented gnuplot-ish text file.
#
# Format version:
# 1.0 -- original release
#-------------------------------------------------------------------------------
import os, time, gc, mpi

from SpheralCompiledPackages import *

class Spheral1dVizDump:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 baseFileName,
                 baseDirectory = ".",
                 listOfFields = [],
                 listOfFieldLists = [],
                 boundaries = []):

        # Store the file name template.
        self.baseFileName = baseFileName

        # Store the base directory name which will be used to create files
        # under.
        self.baseDirectory = baseDirectory
        if self.baseDirectory[-1] == "":
            self.baseDirectory = "."
        if self.baseDirectory[-1] == "/":
            self.baseDirectory = self.baseDirectory[:-1]

        # The internal list of Fields we're dumping.
        self._nodeLists = set()
        self._scalarFieldGroups = {}
        self._vectorFieldGroups = {}
        self._tensorFieldGroups = {}
        self._symTensorFieldGroups = {}

        # If the user volunteered any Fields to be stored, stash them.
        for field in listOfFields:
            self.addField(field)

        # If the user requested any FieldLists be stored, stash their
        # component fields.
        for fieldList in listOfFieldLists:
            self.addFieldList(fieldList)

        # Store the boundaries.
        self._boundaries = boundaries

        # We had better be one dimensional!
        assert len(self._nodeLists) > 0
        for x in self._nodeLists:
            assert "List1d" in str(x)
        self.dimension = "1d"

        # Version of the format written by this class.
        self.version = "1.0"

        return
    
    #---------------------------------------------------------------------------
    # Dump the set of Fields to a dump file.
    #---------------------------------------------------------------------------
    def dump(self, simulationTime, cycle,
             format = "ascii"):

        # Build the name of the directory we will be stuffing the viz file
        # into.
        outputdir = os.path.join(self.baseDirectory, self.baseFileName)

        # Make sure the output directory exists.
        if mpi.rank == 0:
            if not os.path.exists(outputdir):
                try:
                    os.makedirs(outputdir)
                except:
                    raise ValueError("Cannot create output directory %s" % outputdir)
        mpi.barrier()

        # Now build the output file name, including directory.  Make sure
        # the file does not already exist -- if it does we default to overwriting.
        filename = os.path.join(outputdir,
                                self.baseFileName + "-time=%g-cycle=%i.gnu" % (simulationTime, cycle))
##         if os.path.exists(filename):
##             raise ValueError, "File %s already exists!  Aborting." % filename

        # Write the output.
        self.gnuDump(filename,
                     nodeLists = self._nodeLists,
                     time = simulationTime,
                     cycle = cycle,
                     scalarFieldGroups = self._scalarFieldGroups,
                     vectorFieldGroups = self._vectorFieldGroups,
                     tensorFieldGroups = self._tensorFieldGroups,
                     symTensorFieldGroups = self._symTensorFieldGroups)

        # That's it.
        return

    #---------------------------------------------------------------------------
    # Add the given Field to the list of Fields that will be dumped to file.
    #---------------------------------------------------------------------------
    def addField(self, field):
        def addFieldToGroup(field, fieldGroups):
            if field.name in fieldGroups:
                fieldGroups[field.name].append(field)
            else:
                fieldGroups[field.name] = [field]
        if isinstance(field, ScalarField1d):
            addFieldToGroup(field, self._scalarFieldGroups)
        elif isinstance(field, VectorField1d):
            addFieldToGroup(field, self._vectorFieldGroups)
        elif isinstance(field, TensorField1d):
            addFieldToGroup(field, self._tensorFieldGroups)
        elif isinstance(field, SymTensorField1d):
            addFieldToGroup(field, self._symTensorFieldGroups)
        else:
            raise RuntimeError("What is %s?" % field)
        self._nodeLists.add(field.nodeList())
        return

    #---------------------------------------------------------------------------
    # Add the Fields from a FieldList to the set of dumped Fields.
    #---------------------------------------------------------------------------
    def addFieldList(self, fieldList):
        for field in fieldList:
            self.addField(field)
        return

    #---------------------------------------------------------------------------
    # Check if the given NodeList is already listed as being dumped.
    #---------------------------------------------------------------------------
    def haveNodeList(self, nodeList):
        return nodeList in self._nodeLists

    #---------------------------------------------------------------------------
    # Write the text file.
    #---------------------------------------------------------------------------
    def gnuDump(self, filename, nodeLists, time, cycle,
                scalarFieldGroups = {},
                vectorFieldGroups = {},
                tensorFieldGroups = {},
                symTensorFieldGroups = {}):

        def findFieldForNodeList(nodes, fname, flist, FieldConstructor):
            for f in flist:
                if nodes.haveField(f):
                    return f
            return FieldConstructor(fname, nodes)

        # Build up the values for each position.
        values = []
        nodeListID = 0
        nodeListIDs = []
        for nodes in nodeLists:
            nodeListID += 1
            nodeListIDs.append((nodes.name, nodeListID))
            scalarFields = [findFieldForNodeList(nodes, fname, scalarFieldGroups[fname], ScalarField1d) for fname in scalarFieldGroups]
            vectorFields = [findFieldForNodeList(nodes, fname, vectorFieldGroups[fname], VectorField1d) for fname in vectorFieldGroups]
            tensorFields = [findFieldForNodeList(nodes, fname, tensorFieldGroups[fname], TensorField1d) for fname in tensorFieldGroups]
            symTensorFields = [findFieldForNodeList(nodes, fname, symTensorFieldGroups[fname], SymTensorField1d) for fname in symTensorFieldGroups]
            pos = nodes.positions()
            for i in range(nodes.numInternalNodes):
                values.append([pos[i].x, nodeListID])
                for f in scalarFields:
                    values[-1].append(f[i])
                for f in vectorFields:
                    values[-1].append(f[i].x)
                for f in (tensorFields + symTensorFields):
                    values[-1].append(f[i].xx)

        # Gather the values to the root.
        values = mpi.allreduce(values, mpi.SUM)

        # From here on just the root processor is working.
        if mpi.rank == 0:

            # Sort the values by position.
            values.sort()

            # Write the preamble.
            f = open(filename, "w")
            self.writePreamble(f, time, cycle, nodeListIDs,
                               (["pos", "inodelist"] + list(scalarFieldGroups) + list(vectorFieldGroups) + list(tensorFieldGroups) + list(symTensorFieldGroups)))

            # Write the data.
            n = len(values[0])
            for stuff in values:
                assert len(stuff) == n
                f.write((n*" %20g" + "\n") % tuple(stuff))

            f.close()

        return
    
    #---------------------------------------------------------------------------
    # Write the preamble to the given file.
    #---------------------------------------------------------------------------
    def writePreamble(self, f, time, cycle, nodeListIDs, labels):
        f.write("# Time = %g\n# Cycle = %i\n#\n# NodeList -> number mappings:\n" % (time, cycle))
        for (name, id) in nodeListIDs:
            f.write("#    %s  :  %i\n" % (name, id))
        f.write("#")
        for lab in labels:
            f.write(' "%20s"' % lab)
        f.write("\n")
        return

#-------------------------------------------------------------------------------
# Dump out all the Fields in a State object.
# You can pass any of the following for stateThingy:
#     Integrator
#     State
#-------------------------------------------------------------------------------
def dumpPhysicsState(stateThingy,
                     baseFileName,
                     baseDirectory = ".",
                     fields = None,
                     fieldLists = None,
                     currentTime = None,
                     currentCycle = None,
                     dumpGhosts = False,
                     dumpDerivatives = False,
                     boundaries = None):

    # What did we get passed?
    if max([isinstance(stateThingy, x) for x in [Integrator1d, Integrator2d, Integrator3d]]):
        integrator = stateThingy
        dataBase = integrator.dataBase
        state = eval("State%id(integrator.dataBase, integrator.physicsPackages())" % integrator.dataBase.nDim)
        derivs = None
        if dumpDerivatives:
            derivs = eval("StateDerivatives%id(integrator.dataBase, integrator.physicsPackages())" % integrator.dataBase.nDim)
        currentTime = integrator.currentTime
        currentCycle = integrator.currentCycle
    elif max([isinstance(stateThingy, x) for x in [State1d, State2d, State3d]]):
        integrator = None
        dataBase = None
        state = stateThingy
        derivs = None
        assert currentTime is not None
        assert currentCycle is not None
    assert state is not None

    # Did the user specify any data to be dumped?
    if not fields:
        fields = []
    if not fieldLists:
        fieldLists = []

    # Build up the list of fields in the state object.
    fields += [x for x in state.allScalarFields()]
    fields += [x for x in state.allVectorFields()]
    fields += [x for x in state.allTensorFields()]
    fields += [x for x in state.allSymTensorFields()]

    # Are we also dumping the derivative fields?
    if not derivs is None:
        fields += [x for x in derivs.allScalarFields()]
        fields += [x for x in derivs.allVectorFields()]
        fields += [x for x in derivs.allTensorFields()]
        fields += [x for x in derivs.allSymTensorFields()]

    # If available, add the work, H inverse and hmin, hmax, & hmin_hmax_ratio by default.
    if dataBase:
        work = dataBase.globalWork
        fieldLists.append(work)
        Hfl = dataBase.fluidHfield
        Hi = dataBase.newGlobalSymTensorFieldList()
        dataBase.fluidHinverse(Hi)
        fieldLists.append(Hi)
        hmin = dataBase.newGlobalScalarFieldList()
        hmax = dataBase.newGlobalScalarFieldList()
        hminhmax = dataBase.newGlobalScalarFieldList()
        for H, fmin, fmax, fratio in zip(Hfl,
                                         hmin,
                                         hmax,
                                         hminhmax):
            fmin.name = "hmin"
            fmax.name = "hmax"
            fratio.name = "hmin_hmax_ratio"
            n = H.nodeList().numInternalNodes
            for i in range(n):
                ev = H[i].eigenValues()
                fmin[i] = 1.0/ev.maxElement()
                fmax[i] = 1.0/ev.minElement()
                fratio[i] = ev.minElement()/ev.maxElement()
        fieldLists.append(hmin)
        fieldLists.append(hmax)
        fieldLists.append(hminhmax)

        # # We also like to dump the moments of the local point distribution.
        # zerothMoment = eval("ScalarFieldList%id(CopyFields)" % dataBase.nDim)
        # firstMoment = eval("VectorFieldList%id(CopyFields)" % dataBase.nDim)
        # W = eval("TableKernel%id(BSplineKernel%id(), 1000)" % (dataBase.nDim, dataBase.nDim))
        # zerothAndFirstNodalMoments(dataBase.nodeLists, W, True, zerothMoment, firstMoment)
        # fieldLists += [zerothMoment, firstMoment]

    # Add a domain decomposition tag (if we're parallel).
    try:
        import mpi
        domains = dataBase.newGlobalScalarFieldList()
        for f in domains:
            f.name = "Domains"
            for i in range(f.nodeList().numInternalNodes):
                f[i] = mpi.rank
        fieldLists.append(domains)
    except:
        pass

    # Now build the visit dumper.
    dumper = Spheral1dVizDump(baseFileName,
                              baseDirectory,
                              fields,
                              fieldLists,
                              boundaries)

    # Dump the sucker.
    dumper.dump(currentTime, currentCycle)

    # Try to clean up to free memory.
    del dumper, fields, fieldLists
    while gc.collect():
        pass

    return
