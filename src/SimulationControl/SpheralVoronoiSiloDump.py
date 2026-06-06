#-------------------------------------------------------------------------------
# Package for outputting Spheral++ data in silo files on automatically generated
# Voronoi meshes.
# This data can be easily visualized using Visit:
# ftp://ftp.llnl.gov/pub/pub/visit
#
# Format version:
# 1.0 -- original release
#-------------------------------------------------------------------------------
import os, gc, mpi
import time as TIME

from SpheralCompiledPackages import *

from polytope import polytope
from siloMeshDump import *

class SpheralVoronoiSiloDump:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 baseFileName,
                 baseDirectory = ".",
                 listOfFields = [],
                 listOfFieldLists = [],
                 boundaries = [],
                 cells = None,
                 splitCells = False):

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
            self.dimension = "2d"
        elif "List3d" in str(self._nodeLists[0]):
            self.dimension = "3d"
        else:
            assert False

        # Optionally we support the user providing a FieldList of cell geometries.  Not a complete topological mesh,
        # but good enough for writing out.
        if not cells is None:
            assert type(cells) == eval("FacetedVolumeFieldList%s" % self.dimension)
        self.cells = cells
        self.splitCells = splitCells

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
        outputdir = self.baseDirectory # os.path.join(self.baseDirectory, self.baseFileName)

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
                                self.baseFileName + "-time=%g-cycle=%i" % (simulationTime, cycle))
##         if os.path.exists(filename):
##             raise ValueError, "File %s already exists!  Aborting." % filename

        # Did the user provide a FieldList of cell geometries already?
        # start = TIME.clock()
        if self.cells:

            # Yep, so we build a disjoint set of cells as a polytope tessellation.
            mesh = eval("polytope.Tessellation%s()" % self.dimension)
            nDim = eval("Vector%s.nDimensions" % self.dimension)

            # Are we splitting into triangles/tets, or writing the native polygons/polyhera?
            if self.splitCells:
                index2zone = []
                for nodeListi in range(len(self.cells)):
                    n = self.cells[nodeListi].numInternalElements
                    for i in range(n):
                        celli = self.cells(nodeListi, i)
                        verts = celli.vertices()
                        noldnodes = len(mesh.nodes)/nDim
                        noldfaces = len(mesh.faces)
                        noldcells = len(mesh.cells)
                        for j in range(len(verts)):
                            for k in range(nDim):
                                mesh.nodes.append(verts[j][k])

                        if nDim == 2:
                            # Build the triangles
                            PCcelli = PolyClipper.Polygon()
                            PolyClipper.convertToPolygon(PCcelli, celli)
                            tris = PolyClipper.splitIntoTriangles(PCcelli, 1e-10)
                            index2zone.append([])
                            mesh.faces.resize(noldfaces + 3*len(tris))
                            mesh.cells.resize(noldcells + len(tris))
                            for k, tri in enumerate(tris):
                                mesh.faces[noldfaces + 3*k + 0].append(noldnodes + tri[0])
                                mesh.faces[noldfaces + 3*k + 0].append(noldnodes + tri[1])
                                mesh.faces[noldfaces + 3*k + 1].append(noldnodes + tri[1])
                                mesh.faces[noldfaces + 3*k + 1].append(noldnodes + tri[2])
                                mesh.faces[noldfaces + 3*k + 2].append(noldnodes + tri[2])
                                mesh.faces[noldfaces + 3*k + 2].append(noldnodes + tri[0])
                                mesh.cells[noldcells + k].append(noldfaces + 3*k)
                                mesh.cells[noldcells + k].append(noldfaces + 3*k + 1)
                                mesh.cells[noldcells + k].append(noldfaces + 3*k + 2)
                                index2zone[-1].append(noldcells + k)

                        else:
                            # Build the tetrahedra
                            assert nDim == 3
                            PCcelli = PolyClipper.Polyhedron()
                            PolyClipper.convertToPolyhedron(PCcelli, celli)
                            tets = PolyClipper.splitIntoTetrahedra(PCcelli, 1e-10)
                            index2zone.append([])
                            mesh.faces.resize(noldfaces + 4*len(tets))
                            mesh.cells.resize(noldcells + len(tets))
                            for k, tet in enumerate(tets):
                                mesh.faces[noldfaces + 4*k + 0].append(noldnodes + tet[0])
                                mesh.faces[noldfaces + 4*k + 0].append(noldnodes + tet[1])
                                mesh.faces[noldfaces + 4*k + 0].append(noldnodes + tet[2])

                                mesh.faces[noldfaces + 4*k + 1].append(noldnodes + tet[1])
                                mesh.faces[noldfaces + 4*k + 1].append(noldnodes + tet[3])
                                mesh.faces[noldfaces + 4*k + 1].append(noldnodes + tet[2])

                                mesh.faces[noldfaces + 4*k + 2].append(noldnodes + tet[3])
                                mesh.faces[noldfaces + 4*k + 2].append(noldnodes + tet[0])
                                mesh.faces[noldfaces + 4*k + 2].append(noldnodes + tet[2])

                                mesh.faces[noldfaces + 4*k + 3].append(noldnodes + tet[0])
                                mesh.faces[noldfaces + 4*k + 3].append(noldnodes + tet[3])
                                mesh.faces[noldfaces + 4*k + 3].append(noldnodes + tet[1])

                                mesh.cells[noldcells + k].append(noldfaces + 4*k)
                                mesh.cells[noldcells + k].append(noldfaces + 4*k + 1)
                                mesh.cells[noldcells + k].append(noldfaces + 4*k + 2)
                                mesh.cells[noldcells + k].append(noldfaces + 4*k + 3)
                                index2zone[-1].append(noldcells + k)

            else:
                index2zone = None
                copy2polytope(self.cells, mesh)
                # for nodeListi in xrange(len(self.cells)):
                #     n = self.cells[nodeListi].numInternalElements
                #     noldcells = len(mesh.cells)
                #     mesh.cells.resize(noldcells + n)
                #     for i in xrange(n):
                #         celli = self.cells(nodeListi, i)
                #         verts = celli.vertices
                #         facets = celli.facets
                #         noldnodes = len(mesh.nodes)/nDim
                #         noldfaces = len(mesh.faces)
                #         mesh.faces.resize(noldfaces + len(facets))
                #         for j in xrange(len(verts)):
                #             for k in xrange(nDim):
                #                 mesh.nodes.append(verts[j][k])
                #         for j in xrange(len(facets)):
                #             mesh.cells[noldcells + i].append(noldfaces + j)
                #             ipoints = facets[j].ipoints
                #             for k in ipoints:
                #                 mesh.faces[noldfaces + j].append(noldnodes + k)

        else:

            # We need to do the full up polytope tessellation.
            # Build the set of generators from our points.
            gens = vector_of_double()
            nDim = eval("Vector%s.nDimensions" % self.dimension)
            xmin = vector_of_double([1e100]*nDim)
            xmax = vector_of_double([-1e100]*nDim)
            for nodes in self._nodeLists:
                pos = nodes.positions()
                for i in range(nodes.numInternalNodes):
                    for j in range(nDim):
                        gens.append(pos[i][j])
                        xmin[j] = min(xmin[j], pos[i][j])
                        xmax[j] = max(xmax[j], pos[i][j])

            # Check the boundaries for any additional points we want to use for the bounding box.
            for bound in self._boundaries:
                try:
                    pb = dynamicCastBoundaryToPlanarBoundary2d(bound)
                    for p in (pb.enterPlane.point, pb.exitPlane.point):
                        for j in range(nDim):
                            xmin[j] = min(xmin[j], p[j])
                            xmax[j] = max(xmax[j], p[j])
                except:
                    pass

            # Globally reduce and puff up a bit.
            for j in range(nDim):
                xmin[j] = mpi.allreduce(xmin[j], mpi.MIN)
                xmax[j] = mpi.allreduce(xmax[j], mpi.MAX)
                delta = 0.01*(xmax[j] - xmin[j])
                xmin[j] -= delta
                xmax[j] += delta

            # Build the PLC.
            plc = polytope.PLC2d()
            plc.facets.resize(4)
            for i in range(4):
                plc.facets[i].resize(2)
                plc.facets[i][0] = i
                plc.facets[i][1] = (i + 1) % 4
            plccoords = vector_of_double(8)
            plccoords[0] = xmin[0]
            plccoords[1] = xmin[1]
            plccoords[2] = xmax[0]
            plccoords[3] = xmin[1]
            plccoords[4] = xmax[0]
            plccoords[5] = xmax[1]
            plccoords[6] = xmin[0]
            plccoords[7] = xmax[1]

            # Blago!
            # f = open("generators_%i_of_%i.txt" % (mpi.rank, mpi.procs), "w")
            # f.write("# generators x    y\n")
            # for i in xrange(len(gens)/2):
            #     f.write("%g    %g\n" % (gens[2*i], gens[2*i+1]))
            # f.write("# PLC coords    x     y\n")
            # for i in xrange(len(plccoords)/2):
            #     f.write("%g    %g\n" % (plccoords[2*i], plccoords[2*i+1]))
            # f.close()
            # Blago!

            # Build the tessellation.
            if self.dimension == "2d":
                mesh = polytope.Tessellation2d()
                if "TriangleTessellator2d" in dir(polytope):
                    serial_tessellator = polytope.TriangleTessellator2d()
                else:
                    assert "BoostTessellator2d" in dir(polytope)
                    serial_tessellator = polytope.BoostTessellator2d()
            else:
                assert self.dimension == "3d"
                raise RuntimeError("Sorry: 3D tessellation silo dumps are not supported yet.")
            if mpi.procs > 1:
                tessellator = eval("polytope.DistributedTessellator%s(serial_tessellator, False, True)" % self.dimension)
            else:
                tessellator = serial_tessellator
            index2zone = tessellator.tessellateDegenerate(gens, plccoords, plc, 1.0e-8, mesh)

        # print "Took %g sec to generate cells" % (TIME.clock() - start)
        # start = TIME.clock()

        # Figure out how many of each type of field we're dumping.
        intFields = [x for x in self._fields if isinstance(x, eval("IntField%s" % self.dimension))]
                     #[x for x in self._fields if isinstance(x, eval("UnsignedField%s" % self.dimension))] +
                     #[x for x in self._fields if isinstance(x, eval("ULLField%s" % self.dimension))])
        scalarFields = [x for x in self._fields if isinstance(x, eval("ScalarField%s" % self.dimension))]
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
            fvals = f.internalValues()
            for i in range(n.numInternalNodes):
                tr[i] = fvals[i].Trace()
                det[i] = fvals[i].Determinant()
                eigen = fvals[i].eigenValues()
                mineigen[i] = eigen.minElement()
                maxeigen[i] = eigen.maxElement()
            scalarFields += [tr, det, mineigen, maxeigen]
        # print "Took %g sec to build output fields" % (TIME.clock() - start)
        # start = TIME.clock()

        # Write the output.
        timeslice = siloMeshDump(filename, mesh,
                                 index2zone = index2zone,
                                 nodeLists = self._nodeLists,
                                 time = simulationTime,
                                 cycle = cycle,
                                 intFields = intFields,
                                 scalarFields = scalarFields,
                                 vectorFields = vectorFields,
                                 tensorFields = tensorFields,
                                 symTensorFields = symTensorFields)

        # print "Took %g sec to calls siloMeshDump" % (TIME.clock() - start)
        # start = TIME.clock()

        # Write the master file listing all the time slices.
        if mpi.rank == 0:
            mastername = os.path.join(self.baseDirectory, self.masterFileName)
            mf = open(mastername, "a")
            mf.write("%s\n" % timeslice)
            mf.close()
        mpi.barrier()

        # That's it.
        del mesh
        while gc.collect():
            pass
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
                     boundaries = None,
                     splitCells = False):
    assert not dumpGhosts

    # What did we get passed?
    if max([isinstance(stateThingy, x) for x in [Integrator1d, Integrator2d, Integrator3d]]):
        integrator = stateThingy
        dataBase = integrator.dataBase
        state = eval("State%id(integrator.dataBase, integrator.physicsPackages())" % integrator.dataBase.nDim)
        for p in integrator.physicsPackages():
            p.registerAdditionalVisualizationState(dataBase, state)
        derivs = None
        if dumpDerivatives:
            derivs = eval("StateDerivatives%id(integrator.dataBase, integrator.physicsPackages())" % integrator.dataBase.nDim)
        currentTime = integrator.currentTime
        currentCycle = integrator.currentCycle
    elif max([isinstance(stateThingy, x) for x in [State1d, State2d, State3d]]):
        integrator = None
        state = stateThingy
        derivs = None
        assert currentTime is not None
        assert currentCycle is not None
        if isinstance(stateThingy, State1d):
            ndim = 1
        elif isinstance(stateThingy, State2d):
            ndim = 2
        else:
            assert isinstance(stateThingy, State3d)
            ndim = 3
        dataBase = eval("DataBase%id()" % ndim)
        assert state.fieldNameRegistered(HydroFieldNames.mass)
        mass = state.scalarFields(HydroFieldNames.mass)
        for nodes in mass.nodeListPtrs():
            dataBase.appendNodeList(nodes)
        dataBase.updateConnectivityMap(False, False, False)

    assert state is not None and dataBase is not None

    if boundaries is None:
        boundaries = eval("vector_of_Boundary%id()" % dataBase.nDim)

    # Make sure the ghost nodes are set for the Voronoi tessellation work
    for nodes in dataBase.nodeLists:
        nodes.numGhostNodes = 0
        nodes.neighbor().updateNodes()
    for bc in boundaries:
        bc.setAllGhostNodes(dataBase)
        bc.finalizeGhostBoundary()
        for nodes in dataBase.nodeLists:
            nodes.neighbor().updateNodes()
    dataBase.updateConnectivityMap(False, False, False)

    # Did the user specify any data to be dumped?
    if not fields:
        fields = []
    if not fieldLists:
        fieldLists = []

    # Build up the list of fields in the state object.
    fields += [x for x in state.allIntFields()]
    # fields += [x for x in state.allUnsignedFields()]
    # fields += [x for x in state.allULLFields()]
    fields += [x for x in state.allScalarFields()]
    fields += [x for x in state.allVectorFields()]
    fields += [x for x in state.allTensorFields()]
    fields += [x for x in state.allSymTensorFields()]

    # Are we also dumping the derivative fields?
    if not derivs is None:
        fields += [x for x in state.allIntFields()]
        # fields += [x for x in state.allUnsignedFields()]
        # fields += [x for x in state.allULLFields()]
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
        zerothMoment = eval("ScalarFieldList%id(CopyFields)" % dataBase.nDim)
        firstMoment = eval("VectorFieldList%id(CopyFields)" % dataBase.nDim)
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

    # Build the Voronoi-like cells.
    FacetedVolume = {2 : Polygon,
                     3 : Polyhedron}[dataBase.nDim]
    if state.fieldNameRegistered(HydroFieldNames.cells):
        assert state.fieldNameRegistered(HydroFieldNames.cellFaceFlags)
        cells = state.facetedVolumeFields(HydroFieldNames.cells)
        cellFaceFlags = state.vector_of_CellFaceFlagFields(HydroFieldNames.cellFaceFlags)
    else:
        bounds = eval("vector_of_FacetedVolume%id()" % dataBase.nDim)
        holes = eval("vector_of_vector_of_FacetedVolume%id()" % dataBase.nDim)
        weight = eval("ScalarFieldList%id()" % dataBase.nDim)                         # No weights
        surfacePoint = dataBase.newGlobalIntFieldList(0, HydroFieldNames.surfacePoint)
        vol = dataBase.newGlobalScalarFieldList(0.0, HydroFieldNames.volume)
        deltaMedian = dataBase.newGlobalVectorFieldList(eval("Vector%id.zero" % dataBase.nDim), "centroidal delta")
        etaVoidPoints = dataBase.newGlobalvector_of_VectorFieldList(eval("vector_of_Vector%id()" % dataBase.nDim), "eta void points")
        cells = dataBase.newGlobalFacetedVolumeFieldList(FacetedVolume(), "cells")
        cellFaceFlags = dataBase.newGlobalvector_of_CellFaceFlagFieldList(vector_of_CellFaceFlag(), "face flags")
        computeVoronoiVolume(dataBase.globalPosition, 
                             dataBase.globalHfield,
                             dataBase.connectivityMap(),
                             dataBase.solidDamage,
                             bounds,
                             holes,
                             boundaries,
                             weight,
                             surfacePoint,
                             vol,
                             deltaMedian,
                             etaVoidPoints,
                             cells,
                             cellFaceFlags)

    # Amalgamate the cell face flags into a single value per cell.  Not the best visualization yet...
    # cellFaceFlagsSum = dataBase.newGlobalIntFieldList(0, HydroFieldNames.cellFaceFlags + "_sum")
    # for k in xrange(len(cellFaceFlagsSum)):
    #     for i in xrange(len(cellFaceFlagsSum[k])):
    #         cellFaceFlagsSum[k][i] = sum([x.nodeListj for x in cellFaceFlags[k][i]] + [0])
    # fieldLists += [surfacePoint, cellFaceFlagsSum]

    # Now build the visit dumper.
    dumper = SpheralVoronoiSiloDump(baseFileName,
                                    baseDirectory,
                                    fields,
                                    fieldLists,
                                    boundaries,
                                    cells = cells,
                                    splitCells = splitCells)

    # Dump the sucker.
    dumper.dump(currentTime, currentCycle)

    # Try to clean up to free memory.
    del dumper, fields, fieldLists
    while gc.collect():
        pass

    return
