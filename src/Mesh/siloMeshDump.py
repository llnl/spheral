#-------------------------------------------------------------------------------
# Dump Spheral meshed data to a set of Silo files obeying the Overlink 
# specification.
#-------------------------------------------------------------------------------
import os, mpi
from Spheral import *
from SpheralCompiledPackages import silo
from polytope import polytope
import time as TIME

#-------------------------------------------------------------------------------
# A few integer traits by name that Overlink wants.
#-------------------------------------------------------------------------------
OverlinkAttrs = {"ATTR_NODAL"        : 0,
                 "ATTR_ZONAL"        : 1,
                 "ATTR_FACE"         : 2,
                 "ATTR_EDGE"         : 3,
                 "ATTR_INTENSIVE"    : 0,
                 "ATTR_EXTENSIVE"    : 1,
                 "ATTR_FIRST_ORDER"  : 0,
                 "ATTR_SECOND_ORDER" : 1,
                 "ATTR_INTEGER"      : 0,
                 "ATTR_FLOAT"        : 1,
                 }

#-------------------------------------------------------------------------------
# writeSiloMesh -- this is the one the user should actually call!
#-------------------------------------------------------------------------------
def siloMeshDump(dirName, mesh,
                 index2zone = None,
                 nodeLists = [],
                 label = "Spheral++ generated mesh",
                 time = 0.0,
                 cycle = 0,
                 intFields = [],
                 scalarFields = [],
                 vectorFields = [],
                 tensorFields = [],
                 symTensorFields = [],
                 nodeArrays = None,
                 zoneArrays = None,
                 faceArrays = None,
                 pretendRZ = False):

    # print "siloMeshDump: scalarFields: ", [x.name for x in scalarFields]
    # print "              vectorFields: ", [x.name for x in vectorFields]
    # print "              tensorFields: ", [x.name for x in tensorFields]
    # print "           symTensorFields: ", [x.name for x in symTensorFields]
    # print "                 intFields: ", [x.name for x in intFields]

    assert (isinstance(mesh, polytope.Tessellation2d) or
            isinstance(mesh, polytope.Tessellation3d))
    if isinstance(mesh, polytope.Tessellation2d):
        nDim = 2
    else:
        nDim = 3

    # We can only pretend this is an RZ mesh if it's 2D.
    assert (not pretendRZ) or (nDim == 2)
    
    # Make sure the requested directory exists!
    if mpi.rank == 0:
        if not os.path.exists(dirName):
            os.makedirs(dirName)
    mpi.barrier()

    # start = TIME.clock()

    # Extract all the fields we're going to write.
    fieldwad = extractFieldComponents(nodeLists, time, cycle,
                                      intFields, scalarFields, vectorFields, tensorFields, symTensorFields)

    # print "  --> %g sec to extractFieldComponents" % (TIME.clock() - start)
    # start = TIME.clock()

    # index2zone is used to map a Node->(set of zones) if necessary, for instance when tetrahedralizing a polyhedral tessellation of the Nodes.
    if index2zone:
        ntot = sum([nodes.numInternalNodes for nodes in nodeLists])
        assert len(index2zone) == ntot
        ntarget = sum([len(x) for x in index2zone])
        assert ntarget == len(mesh.cells)
        for name, desc, ftype, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
            newsubvars = []
            for subname, vals in subvars:
                valcopy = list(vals)
                vals.resize(ntarget)
                for icells, vali in zip(index2zone, valcopy):
                    for i in icells:
                        vals[i] = vali
                assert len(vals) == ntarget
    else:
        index2zone = [[i,] for i in range(len(mesh.cells))]

    # print "  --> %g sec for index2zone" % (TIME.clock() - start)
    # start = TIME.clock()

    # If we're domain 0 we write the master file.
    masterfile = writeMasterMeshSiloFile(dirName, mesh, label, nodeLists, time, cycle, fieldwad)

    # print "  --> %g sec for masterfile" % (TIME.clock() - start)
    # start = TIME.clock()

    # Each domain writes it's domain file.
    writeDomainMeshSiloFile(dirName, mesh, index2zone, label, nodeLists, time, cycle, fieldwad,
                            pretendRZ, nodeArrays, zoneArrays, faceArrays)

    # print "  --> %g sec for domain file" % (TIME.clock() - start)
    # start = TIME.clock()

    # That's it.
    return masterfile

#-------------------------------------------------------------------------------
# Extract the fields we're going to write.  This requires exploding vector and
# tensor fields into their components.
#-------------------------------------------------------------------------------
def extractFieldComponents(nodeLists, time, cycle,
                           intFields, scalarFields, vectorFields, tensorFields, symTensorFields):
    result = extractFields(nodeLists, time, cycle, intFields,
                           metaDataIntField, extractIntField, dummyIntField)
    result += extractFields(nodeLists, time, cycle, scalarFields,
                            metaDataScalarField, extractScalarField, dummyScalarField)
    result += extractFields(nodeLists, time, cycle, vectorFields,
                            metaDataVectorField, extractVectorField, dummyVectorField)
    result += extractFields(nodeLists, time, cycle, tensorFields,
                            metaDataTensorField, extractTensorField, dummyTensorField)
    result += extractFields(nodeLists, time, cycle, symTensorFields,
                            metaDataTensorField, extractTensorField, dummyTensorField)
    return result

#-------------------------------------------------------------------------------
# Write the master file.
#-------------------------------------------------------------------------------
def writeMasterMeshSiloFile(dirName, mesh, label, nodeLists, time, cycle, fieldwad,
                            meshType = silo.DB_UCDMESH):

    # Figure out who has something to write. Silo can't handle empty variables...
    numZonesPerDomain = mpi.allreduce([(mpi.rank, len(mesh.cells))], mpi.SUM)
    numZonesPerDomain.sort()
    numZonesPerDomain = [x[1] for x in numZonesPerDomain]

    # Only processor 0 actually writes the file.
    linkfile = None
    if mpi.rank == 0:

        nullOpts = silo.DBoptlist()
        
        # Create the master file.
        p0, p1 = os.path.split(dirName)
        fileName = os.path.join(dirName, "OvlTop.silo")
        db = silo.DBCreate(fileName, 
                           silo.DB_CLOBBER, silo.DB_LOCAL, label, silo.DB_HDF5)
        
        # Make directories for variables.
        assert silo.DBMkDir(db, "CELLS") == 0
        assert silo.DBMkDir(db, "POINTS") == 0      # HACK

        # Pattern for constructing per domain variables.
        domainNamePatterns = [("%s/domain%i.silo:" % (p1, i)) + "%s" for i in range(mpi.procs) if numZonesPerDomain[i] > 0]
        numDomains = len(domainNamePatterns)
        
        # Write the domain file names and types.
        domainNames = vector_of_string([p % "MESH" for p in domainNamePatterns])
        meshTypes = vector_of_int([meshType]*numDomains)
        optlist = silo.DBoptlist(8)
        assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
        assert silo.DBPutMultimesh(db, "MMESH", domainNames, meshTypes, optlist) == 0
        
        # Also write out the points as a point mesh.
        if nodeLists and mpi.rank == 0:
            domainNames = vector_of_string([p % "PointMESH" for p in domainNamePatterns])
            meshTypes = vector_of_int([silo.DB_POINTMESH]*numDomains)
            optlist = silo.DBoptlist(8)
            assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
            assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
            assert silo.DBPutMultimesh(db, "MPointMESH", domainNames, meshTypes, optlist) == 0
        
        # Extract the material names, and write per material info if any.
        if nodeLists:
        
            # Write material names (MESH)
            matnames = vector_of_string([x.name for x in nodeLists])
            matnos = vector_of_int(list(range(len(matnames))))
            material_names = vector_of_string([p % "MATERIAL" for p in domainNamePatterns])
            assert len(material_names) == numDomains
            assert len(matnames) == len(nodeLists)
            assert len(matnos) == len(nodeLists)
            optlist = silo.DBoptlist(8)
            assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
            assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
            assert optlist.addOption(silo.DBOPT_MMESH_NAME, "MMESH") == 0
            assert optlist.addOption(silo.DBOPT_MATNAMES, silo.DBOPT_NMATNOS, matnames) == 0
            assert optlist.addOption(silo.DBOPT_MATNOS, silo.DBOPT_NMATNOS, matnos) == 0
            assert silo.DBPutMultimat(db, "MMATERIAL", material_names, optlist) == 0
        
            # Write material names (PointMESH)
            material_names = vector_of_string([p % "PointMATERIAL" for p in domainNamePatterns])
            assert len(material_names) == numDomains
            assert len(matnames) == len(nodeLists)
            assert len(matnos) == len(nodeLists)
            optlist = silo.DBoptlist(8)
            assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
            assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
            assert optlist.addOption(silo.DBOPT_MMESH_NAME, "MPointMESH") == 0
            assert optlist.addOption(silo.DBOPT_MATNAMES, silo.DBOPT_NMATNOS, matnames) == 0
            assert optlist.addOption(silo.DBOPT_MATNOS, silo.DBOPT_NMATNOS, matnos) == 0
            assert silo.DBPutMultimat(db, "MPointMATERIAL", material_names, optlist) == 0
        
            # Write the variable descriptions for non-scalar variables (vector and tensors).
            writeDefvars(db, fieldwad)
        
            # Write the attributes for each variable.
            # This is largely an Overlink thing, which only needs the components of types like vectors.
            if len(fieldwad) > 0:
                attrNames = vector_of_string()
                attrs = vector_of_vector_of_int()
                thpt = vector_of_int([OverlinkAttrs["ATTR_ZONAL"],
                                      OverlinkAttrs["ATTR_INTENSIVE"],
                                      OverlinkAttrs["ATTR_SECOND_ORDER"],
                                      0,
                                      OverlinkAttrs["ATTR_FLOAT"]])
                for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
                    for subname, vals in subvars:
                        attrNames.append(subname)
                        attrs.append(thpt)
                assert len(attrNames) == len(attrs)
                silo.DBPutCompoundarray(db, "VAR_ATTRIBUTES", attrNames, attrs, nullOpts)
        
            # Write the variables descriptors.
            ucdTypes = vector_of_int([silo.DB_UCDVAR]*numDomains)
            for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
                domainVarNames = vector_of_string([p % ("CELLS_" + name) for p in domainNamePatterns])
                assert len(domainVarNames) == numDomains
                assert silo.DBPutMultivar(db, "CELLS/" + name, domainVarNames, ucdTypes, optlistMV) == 0
                if desc != None:
                    for subname, vals in subvars:
                        domainVarNames = vector_of_string([p % ("CELLS_" + subname) for p in domainNamePatterns])
                        assert len(domainVarNames) == numDomains
                        assert silo.DBPutMultivar(db, "CELLS/" + subname, domainVarNames, ucdTypes, optlistVar) == 0
        
            # # HACK: repeating for point values -- remove when vardefs work
            # # Write the variables descriptors.
            # ptTypes = [silo.DB_POINTVAR]*numDomains
            # for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
            #     domainVarNames = []
            #     for p in domainNamePatterns:
            #         domainVarNames.append(p % ("POINTS_" + name))
            #     assert len(domainVarNames) == numDomains
            #     assert silo.DBPutMultivar(db, "POINTS/" + name, domainVarNames, ptTypes, optlistMV) == 0
            #     if desc != None:
            #         for subname, vals in subvars:
            #             domainVarNames = []
            #             for p in domainNamePatterns:
            #                 domainVarNames.append(p % ("POINTS_" + subname))
            #             assert len(domainVarNames) == numDomains
            #             assert silo.DBPutMultivar(db, "POINTS/" + subname, domainVarNames, ptTypes, optlistVar) == 0
        
        # Make a convenient symlink for the master file.
        linkfile = p1 + ".silo"
        fullpath = os.path.join(p0, linkfile)
        if os.path.exists(fullpath):
            os.remove(fullpath)
        targetfile = os.path.join(p1, "OvlTop.silo")
        os.symlink(targetfile, os.path.join(p0, linkfile))

        # That's it.
        assert silo.DBClose(db) == 0

    # Everyone gets the link file name.
    linkfile = mpi.bcast(linkfile, root=0)

    mpi.barrier()
    return linkfile

#-------------------------------------------------------------------------------
# Write the domain file.
#-------------------------------------------------------------------------------
def writeDomainMeshSiloFile(dirName, mesh, index2zone, label, nodeLists, time, cycle, fieldwad,
                            pretendRZ, nodeArrays, zoneArrays, faceArrays,
                            meshType = silo.DB_UCDMESH):

    # Is there anything to do?
    numZones = len(mesh.cells)
    if numZones > 0:

        # Create the file.
        fileName = os.path.join(dirName, "domain%i.silo" % mpi.rank)
        db = silo.DBCreate(fileName, 
                           silo.DB_CLOBBER, silo.DB_LOCAL, label, silo.DB_HDF5)
        nullOpts = silo.DBoptlist()

        # Make directories for variables.
        assert silo.DBMkDir(db, "CELLS") == 0
        assert silo.DBMkDir(db, "POINTS") == 0      # HACK

        # Determine our dimensionality
        if isinstance(mesh, polytope.Tessellation2d):
            nDim = 2
        else:
            assert isinstance(mesh, polytope.Tessellation3d)
            nDim = 3

        # Write a Polygonal zone list.
        zonelistName = { 2 : "zonelist",
                         3 : "PHzonelist" }

        # start = TIME.clock()

        if nDim == 2:
        
            # Read out the zone nodes.  We rely on these already being arranged
            # counter-clockwise.
            zoneNodes = mesh.zoneNodes
            assert len(zoneNodes) == numZones
            assert silo.DBPutZonelist2(db, zonelistName[nDim], nDim, zoneNodes, 0, 0,
                                       vector_of_int([silo.DB_ZONETYPE_POLYGON]*numZones),
                                       vector_of_int([len(zn) for zn in zoneNodes]), # shapesize,
                                       vector_of_int([1]*numZones),
                                       nullOpts) == 0
        
        # Write a Polyhedral zone list.
        if nDim == 3:
        
            # Construct the zone-face list.  We use the ones complement of a face ID
            # to indicate that face needs to be reversed in reference to this zone.
            # This is the same convention as polytope, so just copy it.
            assert silo.DBPutPHZonelist(db, zonelistName[nDim], mesh.facesAsInts, mesh.cells, 0, (numZones - 1), nullOpts) == 0
        
        # print "    --> %g sec to write PHzonelist" % (TIME.clock() - start)
        # start = TIME.clock()

        # Construct the mesh node coordinates.
        assert len(mesh.nodes) % nDim == 0
        if nDim == 2:
            coords = vector_of_vector_of_double([mesh.xnodes, mesh.ynodes])
        else:
            coords = vector_of_vector_of_double([mesh.xnodes, mesh.ynodes, mesh.znodes])
        assert len(coords) == nDim
        
        # print "    --> %g sec to compute coords" % (TIME.clock() - start)
        # start = TIME.clock()

        # Write the mesh itself.
        meshOpts = silo.DBoptlist(8)
        assert meshOpts.addOption(silo.DBOPT_CYCLE, cycle) == 0
        assert meshOpts.addOption(silo.DBOPT_DTIME, time) == 0
        assert meshOpts.addOption(silo.DBOPT_COORDSYS, silo.DB_CARTESIAN) == 0
        assert meshOpts.addOption(silo.DBOPT_NSPACE, nDim) == 0
        assert meshOpts.addOption(silo.DBOPT_TV_CONNECTIVITY, 1) == 0
        if nDim == 2:
            if pretendRZ:
                assert meshOpts.addOption(silo.DBOPT_COORDSYS, silo.DB_CYLINDRICAL) == 0
                assert meshOpts.addOption(silo.DBOPT_XLABEL, "z") == 0
                assert meshOpts.addOption(silo.DBOPT_YLABEL, "r") == 0
            assert silo.DBPutUcdmesh(db, "MESH", coords, numZones, zonelistName[nDim], "NULL", meshOpts) == 0
        else:
            assert meshOpts.addOption(silo.DBOPT_PHZONELIST, zonelistName[nDim]) == 0
            assert silo.DBPutUcdmesh(db, "MESH", coords, numZones, "NULL", "NULL", meshOpts) == 0
        
        # print "    --> %g sec to write mesh" % (TIME.clock() - start)
        # start = TIME.clock()

        # Write materials.
        if nodeLists:
            matnos = vector_of_int(list(range(len(nodeLists))))
            matnames = vector_of_string([nodeList.name for nodeList in nodeLists])
            matlist = []
            offset = 0
            for imat, nodeList in enumerate(nodeLists):
                for i in range(nodeList.numInternalNodes):
                    for j in index2zone[offset + i]:
                        matlist.append(imat)
                offset += nodeList.numInternalNodes
            matlist = vector_of_int(matlist)
            assert len(matlist) == numZones
            assert len(matnames) == len(nodeLists)
            assert len(matnos) == len(nodeLists)
            matOpts = silo.DBoptlist(8)
            assert matOpts.addOption(silo.DBOPT_CYCLE, cycle) == 0
            assert matOpts.addOption(silo.DBOPT_DTIME, time) == 0
            assert matOpts.addOption(silo.DBOPT_MATNAMES, silo.DBOPT_NMATNOS, matnames) == 0
            assert silo.DBPutMaterial(db, "MATERIAL", "MESH", matnos, vector_of_int(matlist), vector_of_int(),
                                      vector_of_int(), vector_of_int(), vector_of_int(), vector_of_double(),
                                      matOpts) == 0
            assert silo.DBPutMaterial(db, "PointMATERIAL", "PointMESH", matnos, vector_of_int(matlist), vector_of_int(),
                                      vector_of_int(), vector_of_int(), vector_of_int(), vector_of_double(),
                                      matOpts) == 0
            # print "      --> %g sec to DBPutMaterial" % (TIME.clock() - start)
            # start = TIME.clock()
        
            # Write the variable descriptions for non-scalar variables (vector and tensors).
            #writeDefvars(db, fieldwad)
        
            # Write the field components.
            varOpts = silo.DBoptlist(8)
            assert varOpts.addOption(silo.DBOPT_CYCLE, cycle) == 0
            assert varOpts.addOption(silo.DBOPT_DTIME, time) == 0
            for name, desc, vtype, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
                for subname, vals in subvars:
                    if len(vals) > 0:
                        if type(vals[0]) == float:
                            ctor = vector_of_double
                        else:
                            ctor = vector_of_int
                        assert silo.DBPutUcdvar1(db, "CELLS_" + subname, "MESH", ctor(vals), ctor([]), silo.DB_ZONECENT, varOpts) == 0

            # print "      --> %g sec to DBPutUcdvar1" % (TIME.clock() - start)
            # start = TIME.clock()

            # # HACK: Write the field components on the point mesh as well.  Remove when the vardef version is working.
            # varOpts = silo.DBoptlist(8)
            # assert varOpts.addOption(silo.DBOPT_CYCLE, cycle) == 0
            # assert varOpts.addOption(silo.DBOPT_DTIME, time) == 0
            # for name, desc, vtype, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
            #     for subname, vals in subvars:
            #         if len(vals) > 0:
            #             if type(vals[0]) == float:
            #                 ctor = vector_of_double
            #             else:
            #                 ctor = vector_of_int
            #             assert silo.DBPutPointvar1(db, "POINTS_" + subname, "PointMESH", ctor(vals), varOpts) == 0

            # print "      --> %g sec to DBPutPointvar1" % (TIME.clock() - start)
            # start = TIME.clock()

        # print "    --> %g sec to write materials" % (TIME.clock() - start)
        # start = TIME.clock()

        # Write the set of neighbor domains.
        thpt = vector_of_vector_of_int()
        thpt.append(vector_of_int([len(mesh.neighborDomains)]))
        thpt.append(vector_of_int([x for x in mesh.neighborDomains]))
        elemNames = vector_of_string(["num neighbor domains", "neighbor domains"])
        assert silo.DBPutCompoundarray(db, "DOMAIN_NEIGHBOR_NUMS", elemNames, thpt, nullOpts) == 0
        
        # Write the shared nodes for each neighbor domain.
        sharedNodes = mesh.sharedNodes
        for ineighborDomain in range(len(mesh.neighborDomains)):
            nodes = [list(sharedNodes[ineighborDomain])]
            assert len(nodes) == len(sharedNodes[ineighborDomain])
            assert silo.DBPutCompoundarray(db, "DOMAIN_NEIGHBOR%i" % ineighborDomain,
                                           ["shared_nodes"],
                                           nodes,
                                           nullOpts) == 0
        
        # If requested, write out annotations for the nodes, zones, and faces.
        if (not (nodeArrays is None) or
            not (zoneArrays is None) or
            not (faceArrays is None)):
            names = []
            values = []
            if not (nodeArrays is None):
                for pair in nodeArrays:
                    assert len(pair) == 2
                    if len(pair[1]) > 0:
                        names.append(pair[0] + "_node")
                        values.append([])
                        for i in pair[1]:
                            values[-1].append(i)
                        assert len(values[-1]) == len(pair[1])
            if not (zoneArrays is None):
                for pair in zoneArrays:
                    assert len(pair) == 2
                    if len(pair[1]) > 0:
                        names.append(pair[0] + "_zone")
                        values.append([])
                        for i in pair[1]:
                            values[-1].append(i)
                        assert len(values[-1]) == len(pair[1])
            if not (faceArrays is None):
                for pair in faceArrays:
                    assert len(pair) == 2
                    if len(pair[1]) > 0:
                        names.append(pair[0] + "_face")
                        values.append([])
                        for i in pair[1]:
                            values[-1].append(i)
                        assert len(values[-1]) == len(pair[1])
            assert len(names) == len(values)
            if len(names) > 0:
                assert silo.DBPutCompoundarray(db, "ANNOTATION_INT", names, values, nullOpts) == 0
        
        # Write the point mesh.
        if nodeLists:
            ntot = sum([n.numInternalNodes for n in nodeLists])
            coords = []
            for j in range(nDim):
                coords.append([])
            for nodes in nodeLists:
                pos = nodes.positions().internalValues()
                n = len(pos)
                for j in range(nDim):
                    for i in range(n):
                        coords[j].append(pos[i][j])
            for j in range(nDim):
                assert len(coords[j]) == ntot
            coords = vector_of_vector_of_double([vector_of_double(vals) for vals in coords])

            # Write the Pointmesh.
            meshOpts = silo.DBoptlist(8)
            assert meshOpts.addOption(silo.DBOPT_CYCLE, cycle) == 0
            assert meshOpts.addOption(silo.DBOPT_DTIME, time) == 0
            assert silo.DBPutPointmesh(db, "PointMESH", coords, meshOpts) == 0

        # That's it.
        assert silo.DBClose(db) == 0

    mpi.barrier()
    return

#-------------------------------------------------------------------------------
# Generic master method to extract full sets of field values.
#-------------------------------------------------------------------------------
def extractFields(nodeLists, time, cycle, fields,
                  metaDataMethod,
                  extractMethod,
                  dudMethod):
    result = []
    if len(fields) > 0:
        dim = dimension(fields[0])

        # Figure out how many values we should have for each field.
        nvals = sum([nodes.numInternalNodes for nodes in nodeLists])

        # Group the fields by name (across all NodeLists).
        fieldsByName = {}
        for f in fields:
            if f.name in fieldsByName:
                fieldsByName[f.name][f.nodeList().name] = f
            else:
                fieldsByName[f.name] = {f.nodeList().name : f}

        # Go through and build up the full set of values.
        for fieldName in fieldsByName:
            # start = TIME.clock()
            subfields = fieldsByName[fieldName]
            assert len(subfields) <= len(nodeLists)

            # Build the meta-data for this field entry.
            name = str(fieldName).replace(" ", "_").replace("-", "_")
            varDef, varType, optlistDef, optlistMV, optlistVar = metaDataMethod(name, time, cycle, dim)

            # Build the complete values for this field.
            vals = []
            for nodes in nodeLists:
                if nodes.name in subfields:
                    vals = extractMethod(name, subfields[nodes.name].internalValues(), vals, dim)
                else:
                    vals = dudMethod(name, nodes.numInternalNodes, vals, dim)
            for thpt in vals:
                assert len(thpt) == 2
                if len(thpt[1]) != nvals:
                    print("BLAGO : ", name, len(thpt[1]), nvals, (nodes.name in subfields))
                assert len(thpt[1]) == nvals

            result.append((name, varDef, varType, optlistDef, optlistMV, optlistVar, vals))
            # print "        **> %g secs for field %s" % (TIME.clock() - start, name)

    return result

#-------------------------------------------------------------------------------
# Int field components.
#-------------------------------------------------------------------------------
def extractIntField(name, field, vals, dim):
    if vals == []:
        vals = [[name, field]]
    else:
        vals[0][1] += field
    return vals

def dummyIntField(name, n, vals, dim):
    if vals == []:
        vals = [[name, [0]*n]]
    else:
        vals[0][1].extend([0]*n)
    return vals

def metaDataIntField(name, time, cycle, dim):
    optlistDef = None
    optlistMV = silo.DBoptlist()
    optlistVar = silo.DBoptlist()
    for optlist in (optlistMV, optlistVar):
        assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
        assert optlist.addOption(silo.DBOPT_TENSOR_RANK, silo.DB_VARTYPE_SCALAR) == 0
    return (None, silo.DB_VARTYPE_SCALAR, optlistDef, optlistMV, optlistVar)

#-------------------------------------------------------------------------------
# Scalar field components.
#-------------------------------------------------------------------------------
def extractScalarField(name, field, vals, dim):
    if vals == []:
        vals = [[name, field]]
    else:
        vals[0][1] += field
    return vals

def dummyScalarField(name, n, vals, dim):
    if vals == []:
        vals = [[name, [0.0]*n]]
    else:
        vals[0][1].extend([0.0]*n)
    return vals

def metaDataScalarField(name, time, cycle, dim):
    optlistDef = None
    optlistMV = silo.DBoptlist()
    optlistVar = silo.DBoptlist()
    for optlist in (optlistMV, optlistVar):
        assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
        assert optlist.addOption(silo.DBOPT_TENSOR_RANK, silo.DB_VARTYPE_SCALAR) == 0
    return (None, silo.DB_VARTYPE_SCALAR, optlistDef, optlistMV, optlistVar)

#-------------------------------------------------------------------------------
# Vector field components.
#-------------------------------------------------------------------------------
def extractVectorField(name, field, vals, dim):
    assert len(vals) == dim or len(vals) == 0
    assert dim in (2,3)

    if dim == 2:
        if vals == []:
            vals = [["%s_x" % name, []],
                    ["%s_y" % name, []]]
        vals[0][1] += [v.x for v in field]
        vals[1][1] += [v.y for v in field]

    else:
        if vals == []:
            vals = [["%s_x" % name, []],
                    ["%s_y" % name, []],
                    ["%s_z" % name, []]]
        vals[0][1] += [v.x for v in field]
        vals[1][1] += [v.y for v in field]
        vals[2][1] += [v.z for v in field]

    return vals

def dummyVectorField(name, n, vals, dim):
    assert len(vals) == dim or len(vals) == 0
    assert dim in (2,3)
    if vals == []:
        if dim == 2:
            vals = [["%s_x" % name, [0.0]*n],
                    ["%s_y" % name, [0.0]*n]]
        else:
            vals = [["%s_x" % name, [0.0]*n],
                    ["%s_y" % name, [0.0]*n],
                    ["%s_z" % name, [0.0]*n]]
    else:
        for i in range(dim):
            vals[i][1].extend([0.0]*n)
    return vals

def metaDataVectorField(name, time, cycle, dim):
    assert dim in (2,3)
    optlistDef = silo.DBoptlist()
    optlistMV = silo.DBoptlist()
    optlistVar = silo.DBoptlist()
    for optlist in (optlistDef, optlistMV, optlistVar):
        assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
    assert optlistMV.addOption(silo.DBOPT_TENSOR_RANK, silo.DB_VARTYPE_VECTOR) == 0
    assert optlistVar.addOption(silo.DBOPT_HIDE_FROM_GUI, 1) == 0
    assert optlistVar.addOption(silo.DBOPT_TENSOR_RANK, silo.DB_VARTYPE_SCALAR) == 0

    if dim == 2:
        return ("{<CELLS/%s_x>, <CELLS/%s_y>}" % (name, name), silo.DB_VARTYPE_VECTOR,
                optlistDef, optlistMV, optlistVar)
    else:
        return ("{<CELLS/%s_x>, <CELLS/%s_y>, <CELLS/%s_z>}" % (name, name, name), silo.DB_VARTYPE_VECTOR,
                optlistDef, optlistMV, optlistVar)

#-------------------------------------------------------------------------------
# Tensor fields components.
#-------------------------------------------------------------------------------
def extractTensorField(name, field, vals, dim):
    assert len(vals) == dim*dim or len(vals) == 0
    assert dim in (2,3)
    if dim == 2:
        if vals == []:
            vals = [["%s_xx" % name, []],
                    ["%s_xy" % name, []],
                    ["%s_yx" % name, []],
                    ["%s_yy" % name, []]]
        vals[0][1] += [t.xx for t in field]
        vals[1][1] += [t.xy for t in field]
        vals[2][1] += [t.yx for t in field]
        vals[3][1] += [t.yy for t in field]
            
    else:
        if vals == []:
            vals = [["%s_xx" % name, []],
                    ["%s_xy" % name, []],
                    ["%s_xz" % name, []],
                    ["%s_yx" % name, []],
                    ["%s_yy" % name, []],
                    ["%s_yz" % name, []],
                    ["%s_zx" % name, []],
                    ["%s_zy" % name, []],
                    ["%s_zz" % name, []]]

        vals[0][1] += [t.xx for t in field]
        vals[1][1] += [t.xy for t in field]
        vals[2][1] += [t.xz for t in field]
        vals[3][1] += [t.yx for t in field]
        vals[4][1] += [t.yy for t in field]
        vals[5][1] += [t.yz for t in field]
        vals[6][1] += [t.zx for t in field]
        vals[7][1] += [t.zy for t in field]
        vals[8][1] += [t.zz for t in field]

    return vals

def dummyTensorField(name, n, vals, dim):
    assert len(vals) == dim*dim or len(vals) == 0
    assert dim in (2,3)
    if vals == []:
        if dim == 2:
            vals = [["%s_xx" % name, [0.0]*n],
                    ["%s_xy" % name, [0.0]*n],
                    ["%s_yx" % name, [0.0]*n],
                    ["%s_yy" % name, [0.0]*n]]
        else:
            vals = [["%s_xx" % name, [0.0]*n],
                    ["%s_xy" % name, [0.0]*n],
                    ["%s_xz" % name, [0.0]*n],
                    ["%s_yx" % name, [0.0]*n],
                    ["%s_yy" % name, [0.0]*n],
                    ["%s_yz" % name, [0.0]*n],
                    ["%s_zx" % name, [0.0]*n],
                    ["%s_zy" % name, [0.0]*n],
                    ["%s_zz" % name, [0.0]*n]]
    else:
        for i in range(dim*dim):
            vals[i][1].extend([0.0]*n)
    return vals

def metaDataTensorField(name, time, cycle, dim):
    assert dim in (2,3)
    optlistDef = silo.DBoptlist()
    optlistMV = silo.DBoptlist()
    optlistVar = silo.DBoptlist()
    for optlist in (optlistDef, optlistMV, optlistVar):
        assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
    assert optlistMV.addOption(silo.DBOPT_TENSOR_RANK, silo.DB_VARTYPE_TENSOR) == 0
    assert optlistVar.addOption(silo.DBOPT_HIDE_FROM_GUI, 1) == 0
    assert optlistVar.addOption(silo.DBOPT_TENSOR_RANK, silo.DB_VARTYPE_SCALAR) == 0

    if dim == 2:
        return ("{{<CELLS/%s_xx>, <CELLS/%s_xy>}, {<CELLS/%s_yx>, <CELLS/%s_yy>}}" % (name, name, name, name), silo.DB_VARTYPE_TENSOR,
                optlistDef, optlistMV, optlistVar)
    else:
        return ("{{<CELLS/%s_xx>, <CELLS/%s_xy>, <CELLS/%s_xz>}, {<CELLS/%s_yx>, <CELLS/%s_yy>, <CELLS/%s_yz>}, {<CELLS/%s_zx>, <CELLS/%s_zy>, <CELLS/%s_zz>}}" % (name, name, name,
                                                                                                                                                                   name, name, name,
                                                                                                                                                                   name, name, name),
                silo.DB_VARTYPE_TENSOR, optlistDef, optlistMV, optlistVar)

#-------------------------------------------------------------------------------
# Write out the prescriptions for defined/derived variables.
#-------------------------------------------------------------------------------
def writeDefvars(db, fieldwad):

    # Write the variable descriptors for non-scalar types (vector and tensor).
    names, defs, types, opts = vector_of_string(), vector_of_string(), vector_of_int(), [] # vector_of_DBoptlist()
    for name, desc, vtype, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
        if desc != None:
            assert optlistDef != None
            assert len(subvars) > 1
            names.append("CELLS/" + name)
            defs.append(desc)
            types.append(vtype)
            opts.append(optlistDef)

    # Map all variables from the MMESH -> MPointMesh.
    showOptlist = silo.DBoptlist()
    hideOptlist = silo.DBoptlist()
    assert showOptlist.addOption(silo.DBOPT_HIDE_FROM_GUI, 0) == 0
    assert hideOptlist.addOption(silo.DBOPT_HIDE_FROM_GUI, 1) == 0
    for name, desc, vtype, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
        if desc == None:
            # This is a simple scalar type
            names.append("POINTS/" + name)
            defs.append('recenter(conn_cmfe(<[0]id:CELLS/%s>,<MPointMESH>))' % name)
            types.append(vtype)
            #types.append(silo.DB_POINTVAR)
            opts.append(showOptlist)

        else:
            # This is a compound type (has components like {x,y,z}
            assert optlistDef != None
            assert len(subvars) > 1
            names.append("POINTS/" + name)
            defs.append(desc.replace("CELLS", "POINTS"))
            types.append(vtype)
            opts.append(optlistDef)

            # Now recenter the scalar components of this compound type
            for subvar in subvars:
                subname = subvar[0]
                names.append("POINTS/" + subname)
                defs.append('recenter(conn_cmfe(<[0]id:CELLS/%s>,<MPointMESH>))' % subname)
                types.append(silo.DB_VARTYPE_SCALAR)
                opts.append(hideOptlist)

    if len(names) > 0:
        assert silo.DBPutDefvars(db, "VARDEFS", names, types, defs, opts) == 0

#-------------------------------------------------------------------------------
# Extract the dimensionality of a field.
#-------------------------------------------------------------------------------
def dimension(field):
    if "Field2d" in str(field):
        return 2
    elif "Field3d" in str(field):
        return 3
    else:
        assert False
