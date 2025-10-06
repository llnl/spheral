#-------------------------------------------------------------------------------
# Dump Spheral point data to a set of Silo files using the Pointmesh silo 
# structures.
#-------------------------------------------------------------------------------
import os, mpi
from Spheral import *
from SpheralCompiledPackages import silo

#-------------------------------------------------------------------------------
# siloPointMeshDump -- this is the one the user should actually call!
#-------------------------------------------------------------------------------
def siloPointmeshDump(baseName, 
                      fields = [],
                      fieldLists = [],
                      baseDirectory = ".",
                      procDirBaseName = "proc-%06i",
                      label = "Spheral++ point mesh",
                      time = 0.0,
                      cycle = 0,
                      dumpGhosts = False):

    # You have to give us something!
    if len(fields) + len(fieldLists) == 0:
        raise ValueError("siloPointmeshDump called with no information to write.")

    # Get the set of NodeLists we're working with.
    nodeListsDict = {}
    for field in fields:
        nodes = field.nodeList()
        nodeListsDict[nodes.name] = nodes
    for fl in fieldLists:
        nls = fl.nodeListPtrs()
        for nodes in nls:
            nodeListsDict[nodes.name] = nodes
    nodeLists = [nodeListsDict[name] for name in nodeListsDict]
    assert len(nodeLists) > 0

    # If needed, create the subdirectories for the processor files.
    if mpi.rank == 0:
        if not os.path.exists(baseDirectory):
            os.makedirs(baseDirectory)
        for i in range(mpi.procs):
            dire = os.path.join(baseDirectory, procDirBaseName % i)
            if not os.path.exists(dire):
                os.makedirs(dire)
    mpi.barrier()

    # We can only pretend this is an RZ mesh if it's 2D.
    ndim = dimension(nodeLists[0])
    if not ndim in (2, 3):
        raise ValueError("You need to provide 2D or 3D information for siloPointMeshDump.")
    
    # Characterize the fields we're going to write.
    allfields = fields[:]
    for fl in fieldLists:
        for f in fl:
            allfields.append(f)
    intFields, scalarFields, vectorFields, tensorFields, symTensorFields = [], [], [], [], []
    for f in allfields:
        if isinstance(f, eval("IntField%id" % ndim)):
            intFields.append(f)
        elif isinstance(f, eval("ScalarField%id" % ndim)):
            scalarFields.append(f)
        elif isinstance(f, eval("VectorField%id" % ndim)):
            vectorFields.append(f)
        elif isinstance(f, eval("TensorField%id" % ndim)):
            tensorFields.append(f)
        elif isinstance(f, eval("SymTensorField%id" % ndim)):
            symTensorFields.append(f)
        else:
            print("siloPointmeshDump WARNING: ignoring unknown field type.")

    # For any tensor fields, dump the trace, determinant, min, and max eigen values.
    for f in (tensorFields + symTensorFields):
        n = f.nodeList()
        tr = eval("ScalarField%id('%s_trace', n)" % (ndim, f.name))
        det = eval("ScalarField%id('%s_determinant', n)" % (ndim, f.name))
        mineigen = eval("ScalarField%id('%s_eigen_min', n)" % (ndim, f.name))
        maxeigen = eval("ScalarField%id('%s_eigen_max', n)" % (ndim, f.name))
        if dumpGhosts:
            nvals = n.numNodes
        else:
            nvals = n.numInternalNodes
        for i in range(nvals):
            eigen = f[i].eigenValues()
            tr[i] = f[i].Trace()
            det[i] = f[i].Determinant()
            mineigen[i] = eigen.minElement()
            maxeigen[i] = eigen.maxElement()
        scalarFields += [tr, det, mineigen, maxeigen]

    # Extract all the fields we're going to write.
    fieldwad = extractFieldComponents(nodeLists, time, cycle, dumpGhosts,
                                      intFields, scalarFields, vectorFields, tensorFields, symTensorFields)

    # If we're domain 0 we write the master file.
    writeMasterSiloFile(ndim, baseDirectory, baseName, procDirBaseName, nodeLists, label, time, cycle, dumpGhosts, fieldwad)

    # Each domain writes it's domain file.
    writeDomainSiloFile(ndim, baseDirectory, baseName, procDirBaseName, nodeLists, label, time, cycle, dumpGhosts, fieldwad)

#-------------------------------------------------------------------------------
# Extract the fields we're going to write.  This requires exploding vector and
# tensor fields into their components.
#-------------------------------------------------------------------------------
def extractFieldComponents(nodeLists, time, cycle, dumpGhosts,
                           intFields, scalarFields, vectorFields, tensorFields, symTensorFields):
    result = extractFields(nodeLists, time, cycle, dumpGhosts, intFields, 
                           metaDataIntField, extractIntField, dummyIntField)
    result += extractFields(nodeLists, time, cycle, dumpGhosts, scalarFields, 
                           metaDataScalarField, extractScalarField, dummyScalarField)
    result += extractFields(nodeLists, time, cycle, dumpGhosts, vectorFields, 
                            metaDataVectorField, extractVectorField, dummyVectorField)
    result += extractFields(nodeLists, time, cycle, dumpGhosts, tensorFields,
                            metaDataTensorField, extractTensorField, dummyTensorField)
    result += extractFields(nodeLists, time, cycle, dumpGhosts, symTensorFields,
                            metaDataTensorField, extractTensorField, dummyTensorField)
    return result

#-------------------------------------------------------------------------------
# Write the master file.
#-------------------------------------------------------------------------------
def writeMasterSiloFile(ndim, baseDirectory, baseName, procDirBaseName, nodeLists,
                        label, time, cycle, dumpGhosts, fieldwad):

    nullOpts = silo.DBoptlist()

    # Make sure which domains have information
    ntot = []
    if dumpGhosts:
        nloc = sum([n.numNodes for n in nodeLists])
    else:
        nloc = sum([n.numInternalNodes for n in nodeLists])
    for iproc in range(mpi.procs):
        ntot.append(mpi.bcast(nloc, iproc))

    # Pattern for constructing per domain variables.
    domainNamePatterns = [os.path.join(procDirBaseName % i, baseName + ".silo:%s") for i in range(mpi.procs) if ntot[i] > 0]
    ndoms = len(domainNamePatterns)

    # Create the master file.
    if mpi.rank == 0:
        fileName = os.path.join(baseDirectory, baseName + ".silo")
        db = silo.DBCreate(fileName, 
                           silo.DB_CLOBBER, silo.DB_LOCAL, label, silo.DB_HDF5)

        # Write the domain file names and types.
        domainNames = vector_of_string()
        meshTypes = vector_of_int([silo.DB_POINTMESH]*ndoms)
        for p in domainNamePatterns:
            domainNames.append(p % "mesh")
        optlist = silo.DBoptlist(1024)
        assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
        assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
        assert silo.DBPutMultimesh(db, "MESH", domainNames, meshTypes, optlist) == 0

    # Extract the material names, and write per material info if any.
    if len(nodeLists) > 0:

        # Write material names.
        if mpi.rank == 0:
            materialNames = [x.name for x in nodeLists]
            material_names = []
            matnames = []
            matnos = []
            for p in domainNamePatterns:
                material_names.append(p % "material")
            for (name, i) in zip(materialNames, list(range(len(materialNames)))):
                matnames.append(name)
                matnos.append(i)
            assert len(material_names) == ndoms
            assert len(matnames) == len(nodeLists)
            assert len(matnos) == len(nodeLists)
            optlist = silo.DBoptlist(1024)
            assert optlist.addOption(silo.DBOPT_CYCLE, cycle) == 0
            assert optlist.addOption(silo.DBOPT_DTIME, time) == 0
            assert optlist.addOption(silo.DBOPT_MATNAMES, silo.DBOPT_NMATNOS, vector_of_string(matnames)) == 0
            assert optlist.addOption(silo.DBOPT_MATNOS, silo.DBOPT_NMATNOS, vector_of_int(matnos)) == 0
            assert silo.DBPutMultimat(db, "MATERIAL", vector_of_string(material_names), optlist) == 0
        
            # Write the variable descriptions for non-scalar variables (vector and tensors).
            writeDefvars(db, fieldwad)

        # Write the variables descriptors.
        types = vector_of_int([silo.DB_POINTVAR]*ndoms)
        for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
            domainVarNames = vector_of_string()
            nlocalvals = sum([len(x[1]) for x in subvars])
            if desc is None:
                for iproc, p in enumerate(domainNamePatterns):
                    nvals = mpi.bcast(nlocalvals, root=iproc)
                    if nvals > 0:
                        domainVarNames.append(p % name)
                    else:
                        domainVarNames.append("EMPTY")
                assert len(domainVarNames) == ndoms
                if mpi.rank == 0:
                    assert silo.DBPutMultivar(db, name, domainVarNames, types, optlistMV) == 0
            else:
                for subname, vals in subvars:
                    domainVarNames = vector_of_string()
                    for iproc, p in enumerate(domainNamePatterns):
                        nvals = mpi.bcast(nlocalvals, root=iproc)
                        if nvals > 0:
                            domainVarNames.append(p % subname)
                        else:
                            domainVarNames.append("EMPTY")
                    assert len(domainVarNames) == ndoms
                    if mpi.rank == 0:
                        assert silo.DBPutMultivar(db, subname, domainVarNames, types, optlistVar) == 0

    # That's it.
    if mpi.rank == 0:
        assert silo.DBClose(db) == 0
        del db

    return

#-------------------------------------------------------------------------------
# Write the domain file.
#-------------------------------------------------------------------------------
def writeDomainSiloFile(ndim, baseDirectory, baseName, procDirBaseName, nodeLists,
                        label, time, cycle, dumpGhosts, fieldwad):

    # How many points are we dumping?
    if dumpGhosts:
        ntot = sum([n.numNodes for n in nodeLists])
    else:
        ntot = sum([n.numInternalNodes for n in nodeLists])

    if ntot:

        # Create the file.
        fileName = os.path.join(baseDirectory, procDirBaseName % mpi.rank, baseName + ".silo")
        db = silo.DBCreate(fileName, 
                           silo.DB_CLOBBER, silo.DB_LOCAL, label, silo.DB_HDF5)
        nullOpts = silo.DBoptlist()

        # Read the per material info.
        coords = [[] for i in range(ndim)]
        for nodes in nodeLists:
            if dumpGhosts:
                pos = nodes.positions().allValues()
            else:
                pos = nodes.positions().internalValues()
            n = len(pos)
            for j in range(ndim):
                for i in range(n):
                    coords[j].append(pos[i][j])
        for j in range(ndim):
            assert len(coords[j]) == ntot
        coords = vector_of_vector_of_double([vector_of_double(icoords) for icoords in coords])

        # Write the Pointmesh.
        meshOpts = silo.DBoptlist(1024)
        assert meshOpts.addOption(silo.DBOPT_CYCLE, cycle) == 0
        assert meshOpts.addOption(silo.DBOPT_DTIME, time) == 0
        assert silo.DBPutPointmesh(db, "mesh", coords, meshOpts) == 0

        # Write materials.
        matnos = [i for i in range(len(nodeLists))]
        assert len(matnos) == len(nodeLists)
        matlist = []
        matnames = []
        for (nodeList, imat) in zip(nodeLists, range(len(nodeLists))):
            if dumpGhosts:
                matlist += [imat]*nodeList.numNodes
            else:
                matlist += [imat]*nodeList.numInternalNodes
            matnames.append(nodeList.name)
        matlist = vector_of_int(matlist)
        matnames = vector_of_string(matnames)
        assert len(matlist) == ntot
        assert len(matnames) == len(nodeLists)
        matOpts = silo.DBoptlist(1024)
        assert matOpts.addOption(silo.DBOPT_CYCLE, cycle) == 0
        assert matOpts.addOption(silo.DBOPT_DTIME, time) == 0
        assert matOpts.addOption(silo.DBOPT_MATNAMES, silo.DBOPT_NMATNOS, vector_of_string(matnames)) == 0
        vecInt = vector_of_int()
        vecDouble = vector_of_double()
        assert silo.DBPutMaterial(db, "material", "mesh", vector_of_int(matnos), matlist, vecInt, 
                                  vecInt, vecInt, vecInt, vecDouble,
                                  matOpts) == 0

        # Write the variable descriptions for non-scalar variables (vector and tensors).
        writeDefvars(db, fieldwad)

        # Write the field components.
        varOpts = silo.DBoptlist(1024)
        assert varOpts.addOption(silo.DBOPT_CYCLE, cycle) == 0
        assert varOpts.addOption(silo.DBOPT_DTIME, time) == 0
        for name, desc, dtype, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
            for subname, vals in subvars:
                if len(vals) > 0:
                    if type(vals[0]) == float:
                        ctor = vector_of_double
                    else:
                        ctor = vector_of_int
                    assert silo.DBPutPointvar1(db, subname, "mesh", ctor(vals), varOpts) == 0

        # That's it.
        assert silo.DBClose(db) == 0
        del db

    return

#-------------------------------------------------------------------------------
# Generic master method to extract full sets of field values.
#-------------------------------------------------------------------------------
def extractFields(nodeLists, time, cycle, dumpGhosts, fields,
                  metaDataMethod,
                  extractMethod,
                  dudMethod):
    result = []
    if len(fields) > 0:
        dim = dimension(fields[0])

        # Figure out how many values we should have for each field.
        if dumpGhosts:
            nvals = sum([nodes.numNodes for nodes in nodeLists])
        else:
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
            subfields = fieldsByName[fieldName]
            assert len(subfields) <= len(nodeLists)

            # Build the meta-data for this field entry.
            name = str(fieldName).replace(" ", "_").replace("-", "_")
            varDef, varType, optlistDef, optlistMV, optlistVar = metaDataMethod(name, time, cycle, dim)

            # Build the complete values for this field.
            vals = []
            for nodes in nodeLists:
                if nodes.name in subfields:
                    if dumpGhosts:
                        vals = extractMethod(name, subfields[nodes.name].allValues(), vals, dim)
                    else:
                        vals = extractMethod(name, subfields[nodes.name].internalValues(), vals, dim)
                else:
                    if dumpGhosts:
                        vals = dudMethod(name, nodes.numNodes, vals, dim)
                    else:
                        vals = dudMethod(name, nodes.numInternalNodes, vals, dim)
            for thpt in vals:
                assert len(thpt) == 2
                assert len(thpt[1]) == nvals

            result.append((name, varDef, varType, optlistDef, optlistMV, optlistVar, vals))

    return result

#-------------------------------------------------------------------------------
# Int field components.
#-------------------------------------------------------------------------------
def extractIntField(name, field, vals, dim):
    if vals == []:
        vals = [[name, field]]
    else:
        vals[0][1].extend(field)
    return vals

def dummyIntField(name, n, vals, dim):
    if vals == []:
        vals = [[name, vector_of_int([0]*n)]]
    else:
        vals[0][1].extend(vector_of_int([0]*n))
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
        vals[0][1].extend(field)
    return vals

def dummyScalarField(name, n, vals, dim):
    if vals == []:
        vals = [[name, vector_of_double([0.0]*n)]]
    else:
        vals[0][1].extend(vector_of_double([0.0]*n))
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
            vals = [["%s_x" % name, vector_of_double()],
                    ["%s_y" % name, vector_of_double()]]
        for v in field:
            vals[0][1].append(v.x)
            vals[1][1].append(v.y)

    else:
        if vals == []:
            vals = [["%s_x" % name, vector_of_double()],
                    ["%s_y" % name, vector_of_double()],
                    ["%s_z" % name, vector_of_double()]]
        for v in field:
            vals[0][1].append(v.x)
            vals[1][1].append(v.y)
            vals[2][1].append(v.z)

    return vals

def dummyVectorField(name, n, vals, dim):
    assert len(vals) == dim or len(vals) == 0
    assert dim in (2,3)
    if vals == []:
        if dim == 2:
            vals = [["%s_x" % name, vector_of_double([0.0]*n)],
                    ["%s_y" % name, vector_of_double([0.0]*n)]]
        else:
            vals = [["%s_x" % name, vector_of_double([0.0]*n)],
                    ["%s_y" % name, vector_of_double([0.0]*n)],
                    ["%s_z" % name, vector_of_double([0.0]*n)]]
    else:
        for i in range(dim):
            vals[i][1].extend(vector_of_double([0.0]*n))
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
        return ("{%s_x, %s_y}" % (name, name), silo.DB_VARTYPE_VECTOR,
                optlistDef, optlistMV, optlistVar)
    else:
        return ("{%s_x, %s_y, %s_z}" % (name, name, name), silo.DB_VARTYPE_VECTOR,
                optlistDef, optlistMV, optlistVar)

#-------------------------------------------------------------------------------
# Tensor fields components.
#-------------------------------------------------------------------------------
def extractTensorField(name, field, vals, dim):
    assert len(vals) == dim*dim or len(vals) == 0
    assert dim in (2,3)
    if dim == 2:
        if vals == []:
            vals = [["%s_xx" % name, vector_of_double()],
                    ["%s_xy" % name, vector_of_double()],
                    ["%s_yx" % name, vector_of_double()],
                    ["%s_yy" % name, vector_of_double()]]
        for t in field:
            vals[0][1].append(t.xx); vals[1][1].append(t.xy)
            vals[2][1].append(t.yx); vals[3][1].append(t.yy)

    else:
        if vals == []:
            vals = [["%s_xx" % name, vector_of_double()],
                    ["%s_xy" % name, vector_of_double()],
                    ["%s_xz" % name, vector_of_double()],
                    ["%s_yx" % name, vector_of_double()],
                    ["%s_yy" % name, vector_of_double()],
                    ["%s_yz" % name, vector_of_double()],
                    ["%s_zx" % name, vector_of_double()],
                    ["%s_zy" % name, vector_of_double()],
                    ["%s_zz" % name, vector_of_double()]]
        for t in field:
            vals[0][1].append(t.xx); vals[1][1].append(t.xy); vals[2][1].append(t.xz)
            vals[3][1].append(t.yx); vals[4][1].append(t.yy); vals[5][1].append(t.yz)
            vals[6][1].append(t.zx); vals[7][1].append(t.zy); vals[8][1].append(t.zz)

    return vals

def dummyTensorField(name, n, vals, dim):
    assert len(vals) == dim*dim or len(vals) == 0
    assert dim in (2,3)
    if vals == []:
        if dim == 2:
            vals = [["%s_xx" % name, vector_of_double([0.0]*n)],
                    ["%s_xy" % name, vector_of_double([0.0]*n)],
                    ["%s_yx" % name, vector_of_double([0.0]*n)],
                    ["%s_yy" % name, vector_of_double([0.0]*n)]]
        else:
            vals = [["%s_xx" % name, vector_of_double([0.0]*n)],
                    ["%s_xy" % name, vector_of_double([0.0]*n)],
                    ["%s_xz" % name, vector_of_double([0.0]*n)],
                    ["%s_yx" % name, vector_of_double([0.0]*n)],
                    ["%s_yy" % name, vector_of_double([0.0]*n)],
                    ["%s_yz" % name, vector_of_double([0.0]*n)],
                    ["%s_zx" % name, vector_of_double([0.0]*n)],
                    ["%s_zy" % name, vector_of_double([0.0]*n)],
                    ["%s_zz" % name, vector_of_double([0.0]*n)]]
    else:
        for i in range(dim*dim):
            vals[i][1].extend(vector_of_double([0.0]*n))
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
        return ("{{%s_xx, %s_xy}, {%s_yx, %s_yy}}" % (name, name, name, name), silo.DB_VARTYPE_TENSOR,
                optlistDef, optlistMV, optlistVar)
    else:
        return ("{{%s_xx, %s_xy, %s_xz}, {%s_yx, %s_yy, %s_yz}, {%s_zx, %s_zy, %s_zz}}" % (name, name, name,
                                                                                           name, name, name,
                                                                                           name, name, name),
                silo.DB_VARTYPE_TENSOR, optlistDef, optlistMV, optlistVar)

#-------------------------------------------------------------------------------
# Write the variable descriptors for non-scalar types (vector and tensor).
#-------------------------------------------------------------------------------
def writeDefvars(db, fieldwad):
    names, defs, types, opts = vector_of_string(), vector_of_string(), vector_of_int(), []
    for name, desc, type, optlistDef, optlistMV, optlistVar, subvars in fieldwad:
        if desc != None:
            assert optlistDef != None
            assert len(subvars) > 1
            names.append(name)
            defs.append(desc)
            types.append(type)
            opts.append(optlistDef)
    if len(names) > 0:
        assert silo.DBPutDefvars(db, "VARDEFS", names, types, defs, opts) == 0

#-------------------------------------------------------------------------------
# Extract the dimensionality of a field.
#-------------------------------------------------------------------------------
def dimension(x):
    dim = int(type(x).__name__[-2:-1])
    return dim
