#-------------------------------------------------------------------------------
# Dumper class for outputting Spheral++ format data for visualization with
# Visit.
# ftp://ftp.llnl.gov/pub/pub/visit
#
# Format version:
# 1.0 -- original release
# 1.1 -- replacing all spaces in strings with underscores, so Visit can
#        manipulate them.
#-------------------------------------------------------------------------------
from SpheralCompiledPackages import *

import os, time, mpi
from siloPointmeshDump import siloPointmeshDump

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
    t0 = time.time()
    dim = type(stateThingy).__name__[-2:]
    if isinstance(stateThingy, eval("Integrator%s" % dim)):
        integrator = stateThingy
        dataBase = integrator.dataBase
        state = eval("State%id(dataBase, integrator.physicsPackages())" % dataBase.nDim)
        for p in integrator.physicsPackages():
            p.registerAdditionalVisualizationState(dataBase, state)
        if dumpDerivatives or dumpGhosts:
            derivs = eval("StateDerivatives%id(dataBase, integrator.physicsPackages())" % dataBase.nDim)
        if dumpGhosts:
            integrator.setGhostNodes()
            integrator.applyGhostBoundaries(state, derivs)
            integrator.finalizeGhostBoundaries()
        currentTime = integrator.currentTime
        currentCycle = integrator.currentCycle
        if not dumpDerivatives:
            derivs = None
    elif isinstance(stateThingy, eval("State%s" % dim)):
        integrator = None
        state = stateThingy
        derivs = None
        assert currentTime is not None
        assert currentCycle is not None
        dataBase = eval("DataBase%s()" % dim)
        assert state.fieldNameRegistered(HydroFieldNames.mass)
        mass = state.scalarFields(HydroFieldNames.mass)
        for nodes in mass.nodeListPtrs():
            dataBase.appendNodeList(nodes)

    assert state is not None
    assert dataBase is not None

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
        fields += [x for x in derivs.allIntFields()]
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
            if dumpGhosts:
                n = H.nodeList().numNodes
            else:
                n = H.nodeList().numInternalNodes
            for i in range(n):
                ev = H[i].eigenValues()
                fmin[i] = 1.0/ev.maxElement()
                fmax[i] = 1.0/max(1e-30, ev.minElement())
                fratio[i] = ev.minElement()/max(1e-30, ev.maxElement())
        fieldLists.append(hmin)
        fieldLists.append(hmax)
        fieldLists.append(hminhmax)

    # Add a domain decomposition tag (if we're parallel).
    try:
        import mpi
        domains = dataBase.newGlobalScalarFieldList()
        for f in domains:
            f.name = "Domains"
            if dumpGhosts:
                n = f.nodeList().numNodes
            else:
                n = f.nodeList().numInternalNodes
            for i in range(n):
                f[i] = mpi.rank
        fieldLists.append(domains)
    except:
        pass

    # Dump the sucker.
    t1 = time.time()
    fullBaseName = baseFileName + "-time=%g-cycle=%i" % (currentTime, currentCycle)
    siloPointmeshDump(baseName = fullBaseName,
                      baseDirectory = baseDirectory,
                      fields = fields,
                      fieldLists = fieldLists,
                      time = currentTime,
                      cycle = currentCycle,
                      dumpGhosts = dumpGhosts)

    # Add to the master file.
    if mpi.rank == 0:
        masterFileName = os.path.join(baseDirectory, baseFileName + ".visit")
        mf = open(masterFileName, "a")
        mf.write("%s\n" % (fullBaseName + ".silo"))
        mf.close()

    t2 = time.time()
    print("SpheralPointMeshSiloDump: spent %g seconds on preliminaries, %g writing files." % (t1 - t0, t2 - t1))

    return

