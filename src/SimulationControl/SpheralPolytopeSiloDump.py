#-------------------------------------------------------------------------------
# Package for outputting Spheral++ data in silo files on automatically generated
# Voronoi meshes.  This version uses polytope internals to do all the heavy 
# lifting.
# This data can be easily visualized using Visit:
# ftp://ftp.llnl.gov/pub/pub/visit
#
# Format version:
# 1.0 -- original release
#-------------------------------------------------------------------------------
import os, sys, time, gc, mpi

from SpheralCompiledPackages import *
from PolytopeModules import polytope

class SpheralPolytopeSiloDump:

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
        self.masterFileName = baseFileName + ".visit"

        # Store the base directory name which will be used to create files
        # under.
        self.baseDirectory = baseDirectory
        if self.baseDirectory[-1] == "":
            self.baseDirectory = "."
        if self.baseDirectory[-1] == "/":
            self.baseDirectory = self.baseDirectory[:-1]

        # The internal list of Fields we're dumping.
        self._fields = []

        # If the user volunteered any Fields to be stored, stash them.
        for field in listOfFields:
            self.addField(field)

        # If the user requested any FieldLists be stored, stash their
        # component fields.
        for fieldList in listOfFieldLists:
            self.addFieldList(fieldList)

        # Build the set of NodeLists for all Fields we're dumping.
        self._nodeLists = self.uniqueNodeLists(self._fields)
        assert len(self._nodeLists) > 0 and len(self._nodeLists) <= len(self._fields)

        # Store the boundaries.
        self._boundaries = boundaries

        # Determine our dimensionality.
        assert len(self._nodeLists) > 0
        if "List2d" in str(self._nodeLists[0]):
            self.ndim = 2
            self.dimension = "2d"
        elif "List3d" in str(self._nodeLists[0]):
            self.ndim = 3
            self.dimension = "3d"
        else:
            assert False

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

        # Now build the output file name.
        filename = self.baseFileName + "-time=%g-cycle=%i" % (simulationTime, cycle)

        # Build the set of generators from our points.
        gens = vector_of_double()
        for nodes in self._nodeLists:
            pos = nodes.positions()
            for i in range(nodes.numInternalNodes):
                for coord in pos[i]:
                    gens.append(coord)

        # Build the tessellation.
        if self.ndim == 2:
            mesh = polytope.Tessellation2d()
            if "TriangleTessellator2d" in dir(polytope):
                serial_tessellator = polytope.TriangleTessellator2d()
            else:
                assert "BoostTessellator2d" in dir(polytope)
                serial_tessellator = polytope.BoostTessellator2d()
        else:
            assert self.ndim == 3
            raise RuntimeError("Sorry: 3D tessellation silo dumps are not supported yet.")
        if mpi.procs > 1:
            tessellator = eval("polytope.DistributedTessellator%s(serial_tessellator)" % self.dimension)
        else:
            tessellator = serial_tessellator
        tessellator.tessellate(gens, mesh)

        # Figure out how many of each type of field we're dumping.
        scalarFields = ([x for x in self._fields if isinstance(x, eval("ScalarField%s" % self.dimension))] +
                        [x for x in self._fields if isinstance(x, eval("IntField%s" % self.dimension))])
        vectorFields = [x for x in self._fields if isinstance(x, eval("VectorField%s" % self.dimension))]
        tensorFields = [x for x in self._fields if isinstance(x, eval("TensorField%s" % self.dimension))]
        symTensorFields = [x for x in self._fields if isinstance(x, eval("SymTensorField%s" % self.dimension))]

        # For tensor fields we like to dump out some extra info.
        for f in (tensorFields + symTensorFields):
            n = f.nodeList()
            tr = eval("ScalarField%s('%s_trace', n)" % (self.dimension, f.name))
            det = eval("ScalarField%s('%s_determinant', n)" % (self.dimension, f.name))
            mineigen = eval("ScalarField%s('%s_eigen_min', n)" % (self.dimension, f.name))
            maxeigen = eval("ScalarField%s('%s_eigen_max', n)" % (self.dimension, f.name))
            for i in range(n.numInternalNodes):
                tr[i] = f[i].Trace()
                det[i] = f[i].Determinant()
                eigen = f[i].eigenValues()
                mineigen[i] = eigen.minElement()
                maxeigen[i] = eigen.maxElement()
            scalarFields += [tr, det, mineigen, maxeigen]

        # Build the dictionary of (name, array) values that the polytope silo method can take.
        # Scalar fields.
        cellFields = {}
        for field in scalarFields:
            cellFields[self.siloMangleName(field.name)] = [float(x) for x in field.internalValues()]

        # Vector fields.
        for field in vectorFields:
            cellFields[field.name + "_x"] = [x.x for x in field.internalValues()]
            cellFields[field.name + "_y"] = [x.y for x in field.internalValues()]
            if self.ndim == 3:
                cellFields[field.name + "_z"] = [x.z for x in field.internalValues()]

        # Tensor fields.
        for field in tensorFields:
            cellFields[field.name + "_xx"] = [x.xx for x in field.internalValues()]
            cellFields[field.name + "_xy"] = [x.xy for x in field.internalValues()]
            cellFields[field.name + "_yx"] = [x.yx for x in field.internalValues()]
            cellFields[field.name + "_yy"] = [x.yy for x in field.internalValues()]
            if self.ndim == 3:
                cellFields[field.name + "_xz"] = [x.xz for x in field.internalValues()]
                cellFields[field.name + "_yz"] = [x.yz for x in field.internalValues()]
                cellFields[field.name + "_zx"] = [x.zx for x in field.internalValues()]
                cellFields[field.name + "_zy"] = [x.zy for x in field.internalValues()]
                cellFields[field.name + "_zz"] = [x.zz for x in field.internalValues()]

        # SymTensor fields.
        for field in symTensorFields:
            cellFields[field.name + "_xx"] = [x.xx for x in field.internalValues()]
            cellFields[field.name + "_xy"] = [x.xy for x in field.internalValues()]
            cellFields[field.name + "_yy"] = [x.yy for x in field.internalValues()]
            if self.ndim == 3:
                cellFields[field.name + "_xz"] = [x.xz for x in field.internalValues()]
                cellFields[field.name + "_yz"] = [x.yz for x in field.internalValues()]
                cellFields[field.name + "_zz"] = [x.zz for x in field.internalValues()]

        # Write the output.
        writeTessellation = eval("polytope.writeTessellation%s" % self.dimension)
        writeTessellation(mesh, filename, outputdir,
                          nodeFields = None,
                          edgeFields = None,
                          faceFields = None,
                          cellFields = cellFields,
                          cycle = cycle,
                          time = simulationTime)

        # Write the master file listing all the time slices.
        timeslice = "%s-%d.silo" % (filename, cycle)
        if mpi.rank == 0:
            mastername = os.path.join(self.baseDirectory, self.baseFileName, self.masterFileName)
            mf = open(mastername, "a")
            mf.write("%s\n" % timeslice)
            mf.close()
        mpi.barrier()

        # That's it.
        return

    #---------------------------------------------------------------------------
    # Add the given Field to the list of Fields that will be dumped to file.
    #---------------------------------------------------------------------------
    def addField(self, field):
        if ((field.name, field.nodeList().name) not in
            [(f.name, f.nodeList().name) for f in self._fields]):
            self._fields.append(field)
        return

    #---------------------------------------------------------------------------
    # Add the Fields from a FieldList to the set of dumped Fields.
    #---------------------------------------------------------------------------
    def addFieldList(self, fieldList):
        for field in fieldList:
            self.addField(field)
        return

    #---------------------------------------------------------------------------
    # Return the unique set of NodeLists based on the Fields we're dumping.
    #---------------------------------------------------------------------------
    def uniqueNodeLists(self, listOfFields):
        redundantSet = [(f.nodeList().name, f.nodeList()) for f in listOfFields]
        redundantSet.sort()
        assert len(redundantSet) > 0
        result = [redundantSet[0][1]]
        for i in range(1, len(redundantSet)):
            if redundantSet[i][0] != redundantSet[i - 1][0]:
                result.append(redundantSet[i][1])
        names = [x.name for x in result]
        assert max([names.count(x) for x in names]) == 1
        return result

    #---------------------------------------------------------------------------
    # Check if the given NodeList is already listed as being dumped.
    #---------------------------------------------------------------------------
    def haveNodeList(self, nodeList):
        return nodeList.name in [n.name for n in self._nodeLists]

    #---------------------------------------------------------------------------
    # Mangle a string into something silo acceptable.
    #---------------------------------------------------------------------------
    def siloMangleName(self, name):
        return name.replace(" ", "_")

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
    assert not dumpGhosts

    # What did we get passed?
    if max([isinstance(stateThingy, x) for x in [Integrator1d, Integrator2d, Integrator3d]]):
        integrator = stateThingy
        dataBase = integrator.dataBase()
        state = eval("State%id(integrator.dataBase(), integrator.physicsPackages())" % integrator.dataBase().nDim)
        derivs = None
        if dumpDerivatives:
            derivs = eval("StateDerivatives%id(integrator.dataBase(), integrator.physicsPackages())" % integrator.dataBase().nDim)
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
    fields += [x for x in state.allIntFields()]
    fields += [x for x in state.allScalarFields()]
    fields += [x for x in state.allVectorFields()]
    fields += [x for x in state.allTensorFields()]
    fields += [x for x in state.allSymTensorFields()]

    # Are we also dumping the derivative fields?
    if not derivs is None:
        fields += [x for x in state.allIntFields()]
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

        # We also like to dump the moments of the local point distribution.
        zerothMoment = eval("ScalarFieldList%id(1)" % dataBase.nDim)
        firstMoment = eval("VectorFieldList%id(1)" % dataBase.nDim)
        W = eval("TableKernel%id(BSplineKernel%id(), 1000)" % (dataBase.nDim, dataBase.nDim))
        zerothAndFirstNodalMoments(dataBase.nodeLists, W, True, zerothMoment, firstMoment)
        fieldLists += [zerothMoment, firstMoment]

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
    dumper = SpheralPolytopeSiloDump(baseFileName,
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
