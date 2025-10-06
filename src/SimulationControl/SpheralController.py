#------------------------------------------------------------------------------
# A simple class to control simulation runs for Spheral.
#------------------------------------------------------------------------------
import sys, os, gc, warnings, mpi

from SpheralCompiledPackages import *
from SpheralTimer import SpheralTimer
from SpheralUtilities import adiak_value, TimerMgr
from SpheralConservation import SpheralConservation
from GzipFileIO import GzipFileIO
from SpheralTestUtilities import globalFrame
from NodeGeneratorBase import ConstantRho
from findLastRestart import findLastRestart

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

class SpheralController:

    #--------------------------------------------------------------------------
    # Constuctor.
    #--------------------------------------------------------------------------
    def __init__(self, integrator,
                 kernel = None,
                 statsStep = 1,
                 printStep = 1,
                 garbageCollectionStep = 100,
                 redistributeStep = None,
                 restartStep = None,
                 restartBaseName = "restart",
                 restartObjects = [],
                 restartFileConstructor = SiloFileIO,
                 SPIOFileCountPerTimeslice = None,
                 restoreCycle = None,
                 initializeDerivatives = False,
                 vizBaseName = None,
                 vizDir =  ".",
                 vizStep = None,
                 vizTime = None,
                 vizMethod = None,
                 vizFields = [],
                 vizFieldLists = [],
                 vizGhosts = False,
                 vizDerivs = False,
                 initialTime = 0.0,
                 SPH = False,
                 periodicWork = [],
                 periodicTimeWork = [],
                 skipInitialPeriodicWork = False,
                 iterateInitialH = True,
                 numHIterationsBetweenCycles = 0,
                 reinitializeNeighborsStep = 10,
                 volumeType = RKVolumeType.RKVoronoiVolume,
                 facetedBoundaries = None,
                 printAllTimers = False):
        self.restartBaseName = restartBaseName
        self.restart = RestartableObject(self)
        self.integrator = integrator
        self.restartObjects = restartObjects
        self.restartFileConstructor = restartFileConstructor
        self.SPIOFileCountPerTimeslice = SPIOFileCountPerTimeslice
        self.numHIterationsBetweenCycles = numHIterationsBetweenCycles
        self._break = False

        # test of dem if so throw in a default kernel (all we need is the extent)
        self.isDEM = self.packagesContainDEM(integrator.physicsPackages())
        if self.isDEM and kernel is None:
            kernel = eval("TableKernel{0}d(WendlandC2Kernel{0}d(), 100)".format(integrator.dataBase.nDim))

        # Extract the interpolation kernel for iterating H and such
        if kernel is None:
            for package in integrator.physicsPackages():
                if hasattr(package, "kernel"):
                    kernel = package.kernel
                    break
        if kernel is None:
            raise RuntimeError("SpheralController: unable to extract an appropriate interpolation kernel, please provide in constructor arguments")
        if isinstance(kernel, SphericalKernel):
            self.kernel = kernel.baseKernel1d
        else:
            self.kernel = kernel

        # Determine the dimensionality of this run, based on the integrator.
        self.dim = "%id" % self.integrator.dataBase.nDim
        adiak_value("dim", self.dim)

        # Determine the visualization method.
        if self.dim == "1d":
            from Spheral1dVizDump import dumpPhysicsState
        else:
            from SpheralVoronoiSiloDump import dumpPhysicsState
        if vizMethod:
            dumpPhysicsState = vizMethod
        self.vizMethod = dumpPhysicsState
        self.vizGhosts = vizGhosts
        self.vizDerivs = vizDerivs

        # Organize the physics packages as appropriate
        self.organizePhysicsPackages(self.kernel, volumeType, facetedBoundaries)

        # If this is a parallel run, automatically construct and insert
        # a DistributedBoundaryCondition into each physics package.
        self.insertDistributedBoundary(integrator.physicsPackages())

        # Should we look for the last restart set?
        if restartBaseName:
            if restoreCycle == -1:
                restoreCycle = findLastRestart(restartBaseName)
        else:
            restoreCycle = None

        # Generic initialization work.
        self.reinitializeProblem(restartBaseName,
                                 vizBaseName,
                                 initialTime = initialTime,
                                 statsStep = statsStep,
                                 printStep = printStep,
                                 garbageCollectionStep = garbageCollectionStep,
                                 redistributeStep = redistributeStep,
                                 restartStep = restartStep,
                                 restoreCycle = restoreCycle,
                                 initializeDerivatives = initializeDerivatives,
                                 vizDir = vizDir,
                                 vizStep = vizStep,
                                 vizTime = vizTime,
                                 vizFields = vizFields,
                                 vizFieldLists = vizFieldLists,
                                 periodicWork = periodicWork,
                                 periodicTimeWork = periodicTimeWork,
                                 skipInitialPeriodicWork = skipInitialPeriodicWork,
                                 iterateInitialH = iterateInitialH,
                                 reinitializeNeighborsStep = reinitializeNeighborsStep,
                                 volumeType = volumeType,
                                 facetedBoundaries = facetedBoundaries)

        # Read the restart information if requested.
        if not restoreCycle is None:
            print("Reading from restart file ", restoreCycle)
            self.loadRestartFile(restoreCycle)
        
        return

    #--------------------------------------------------------------------------
    # Destructor
    #--------------------------------------------------------------------------
    def __del__(self):
        pass

    #--------------------------------------------------------------------------
    # (Re)initialize the current problem (and controller state).
    # This method is intended to be called before the controller begins a new
    # problem from time 0.
    #--------------------------------------------------------------------------
    def reinitializeProblem(self, restartBaseName, vizBaseName,
                            initialTime = 0.0,
                            statsStep = 1,
                            printStep = 1,
                            garbageCollectionStep = 100,
                            redistributeStep = None,
                            restartStep = None,
                            restoreCycle = None,
                            initializeDerivatives = False,
                            vizDir = None,
                            vizStep = None,
                            vizTime = None,
                            vizFields = [],
                            vizFieldLists = [],
                            periodicWork = [],
                            periodicTimeWork = [],
                            skipInitialPeriodicWork = False,
                            iterateInitialH = True,
                            reinitializeNeighborsStep = 10,
                            volumeType = RKVolumeType.RKVoronoiVolume,
                            facetedBoundaries = None):

        # Call the global C++ initialization method
        setGlobalFlags()

        # Intialize the cycle count.
        self.totalSteps = 0

        # Construct a timer to track the cycle step time.
        self.stepTimer = SpheralTimer("Time per integration cycle.")

        # Prepare an empty set of periodic work.
        self._periodicWork = []
        self._periodicTimeWork = []
        
        # Set the restart file base name.
        if restartBaseName:
            self.setRestartBaseName(restartBaseName)
        
        # Set the simulation time.
        self.integrator.currentTime = initialTime

        # Check if we have any boundary conditions that need to copy initial state
        uniquebcs = self.integrator.uniqueBoundaryConditions()
        stateBCs = [eval("InflowOutflowBoundary%id" % i) for i in dims] + [eval("ConstantBoundary%id" % i) for i in dims]
        stateBCactive = max([False] + [isinstance(x, y) for y in stateBCs for x in uniquebcs])
        for bc in uniquebcs:
            bc.initializeProblemStartup(False)

        # Create ghost nodes for the physics packages to initialize with.
        db = self.integrator.dataBase
        db.reinitializeNeighbors()
        self.integrator.setGhostNodes()
        db.updateConnectivityMap(False)

        # Initialize the integrator and packages.
        packages = self.integrator.physicsPackages()
        for package in packages:
            package.initializeProblemStartup(db)
        state = eval("State%s(db, packages)" % (self.dim))
        derivs = eval("StateDerivatives%s(db, packages)" % (self.dim))

        # Build the connectivity
        requireConnectivity = max([pkg.requireConnectivity() for pkg in packages])
        if requireConnectivity:
            requireGhostConnectivity = max([pkg.requireGhostConnectivity() for pkg in packages])
            requireOverlapConnectivity = max([pkg.requireOverlapConnectivity() for pkg in packages])
            requireIntersectionConnectivity = max([pkg.requireIntersectionConnectivity() for pkg in packages])
            db.reinitializeNeighbors()
            db.updateConnectivityMap(requireGhostConnectivity, requireOverlapConnectivity, requireIntersectionConnectivity)
            state.enrollConnectivityMap(db.connectivityMapPtr(requireGhostConnectivity, requireOverlapConnectivity, requireIntersectionConnectivity))

        # Initialize dependent state
        for package in packages:
            package.initializeProblemStartupDependencies(db, state, derivs)

        # If we're starting from scratch, initialize the H tensors.
        if restoreCycle is None and not skipInitialPeriodicWork and iterateInitialH:
            self.iterateIdealH()

        db.reinitializeNeighbors()
        self.integrator.setGhostNodes()
        db.updateConnectivityMap(False)
        self.integrator.applyGhostBoundaries(state, derivs)
        for bc in uniquebcs:
            bc.initializeProblemStartup(False)

        # Reconstruct the state and derivs since H and connectivity have updated
        state = eval("State%s(db, packages)" % (self.dim))
        derivs = eval("StateDerivatives%s(db, packages)" % (self.dim))

        # If requested, initialize the derivatives.
        if initializeDerivatives or stateBCactive:
            self.integrator.preStepInitialize(state, derivs)
            dt = self.integrator.selectDt(self.integrator.dtMin, self.integrator.dtMax, state, derivs)
            self.integrator.initializeDerivatives(initialTime, dt, state, derivs)
            self.integrator.evaluateDerivatives(initialTime, dt, db, state, derivs)

        # If there are stateful boundaries present, give them one more crack at copying inital state
        db.reinitializeNeighbors()
        self.integrator.setGhostNodes()
        db.updateConnectivityMap(False)
        self.integrator.applyGhostBoundaries(state, derivs)
        for bc in uniquebcs:
            bc.initializeProblemStartup(True)
        db.reinitializeNeighbors()
        self.integrator.setGhostNodes()
        db.updateConnectivityMap(False)
        self.integrator.applyGhostBoundaries(state, derivs)

        # Set up the default periodic work.
        self.appendPeriodicWork(self.printCycleStatus, printStep)
        self.appendPeriodicWork(self.garbageCollection, garbageCollectionStep)
        self.appendPeriodicWork(self.updateConservation, statsStep)
        self.appendPeriodicWork(self.updateRestart, restartStep)
        self.appendPeriodicWork(self.reinitializeNeighbors, reinitializeNeighborsStep)
        for x, freq in periodicWork:
            self.appendPeriodicWork(x, freq)
        for x, freq in periodicTimeWork:
            self.appendPeriodicTimeWork(x, freq)

        # Add the dynamic redistribution object to the controller.
        self.addRedistributeNodes(self.kernel)

        # Set up any visualization required.
        if not vizBaseName is None:
            assert not vizDir is None
            assert not (vizStep is None and vizTime is None)
            self.addVisualizationDumps(vizBaseName, vizDir, vizStep, vizTime,
                                       vizFields, vizFieldLists)

        # Construct a fresh conservation check object.
        # Hopefully by this time all packages have initialized their own extra energy bins.
        self.conserve = SpheralConservation(self.integrator.dataBase,
                                            self.integrator.physicsPackages())

        # Force the periodic work to fire at problem initalization.
        if (not skipInitialPeriodicWork) and (restoreCycle is None):
            self.doPeriodicWork(force=True)

        # We add this one after forcing periodic work so it's not always fired right at the beginning of a calculation.
        self.appendPeriodicWork(self.updateDomainDistribution, redistributeStep)

        return

    #--------------------------------------------------------------------------
    # Set the restart base name.
    #--------------------------------------------------------------------------
    def setRestartBaseName(self,
                           name,
                           rank = mpi.rank,
                           procs = mpi.procs):
        self.restartBaseName = name

        # If we're running parallel then add the domain info to the restart
        # base name.
        if procs > 1:
            self.restartBaseName += '_rank%i_of_%idomains' % (rank, procs)

        return

    #--------------------------------------------------------------------------
    # Return the current time.
    #--------------------------------------------------------------------------
    def time(self):
        return self.integrator.currentTime

    #--------------------------------------------------------------------------
    # Return the last timestep.
    #--------------------------------------------------------------------------
    def lastDt(self):
        return self.integrator.lastDt

    #--------------------------------------------------------------------------
    # Smooth the physical variables.
    #--------------------------------------------------------------------------
    def smoothState(self, smoothIters=1):
        db = self.integrator.dataBase
        scalarSmooth = eval("smoothScalarFields%id" % db.nDim)
        vectorSmooth = eval("smoothVectorFields%id" % db.nDim)
        tensorSmooth = eval("smoothSymTensorFields%id" % db.nDim)
        for iter in range(smoothIters):
            state = eval("State%id(db, self.integrator.physicsPackages())" % db.nDim)
            derivs = eval("StateDerivatives%id(db, self.integrator.physicsPackages())" % db.nDim)
            self.integrator.setGhostNodes()
            self.integrator.applyGhostBoundaries(state, derivs)
            smoothedVelocity = vectorSmooth(db.fluidVelocity,
                                            db.fluidPosition,
                                            db.fluidWeight,
                                            db.fluidMass,
                                            db.fluidMassDensity,
                                            db.fluidHfield,
                                            self.kernel)
            smoothedSpecificThermalEnergy = scalarSmooth(db.fluidSpecificThermalEnergy,
                                                         db.fluidPosition,
                                                         db.fluidWeight,
                                                         db.fluidMass,
                                                         db.fluidMassDensity,
                                                         db.fluidHfield,
                                                         self.kernel)
            smoothedHfield = tensorSmooth(db.fluidHfield,
                                          db.fluidPosition,
                                          db.fluidWeight,
                                          db.fluidMass,
                                          db.fluidMassDensity,
                                          db.fluidHfield,
                                          self.kernel)
            db.fluidVelocity.assignFields(smoothedVelocity)
            db.fluidSpecificThermalEnergy.assignFields(smoothedSpecificThermalEnergy)
            db.fluidHfield.assignFields(smoothedHfield)
            for nodeList in db.fluidNodeLists:
                nodeList.neighbor().updateNodes()
            db.updateConnectivityMap()
            for nodeList in db.fluidNodeLists:
                nodeList.updateWeight(db.connectivityMap())

        return

    #--------------------------------------------------------------------------
    # Allow Spheral to smoothly exit out of the controller loop.
    #--------------------------------------------------------------------------
    def stop(self):
        self._break = True
        return

    #--------------------------------------------------------------------------
    # Advance the system to the given simulation time.  The user can also
    # specify a max number of steps to take.
    #--------------------------------------------------------------------------
    def advance(self, goalTime, maxSteps=None):
        TimerMgr.timer_begin("advance")
        currentSteps = 0
        while (self.time() < goalTime and
               (maxSteps == None or currentSteps < maxSteps) and (self._break == False)):
            self.stepTimer.start()
            self.integrator.step(goalTime)
            if self.numHIterationsBetweenCycles > 0:
                self.iterateIdealH(maxIdealHIterations = self.numHIterationsBetweenCycles)
            self.stepTimer.stop()
            currentSteps = currentSteps + 1
            self.totalSteps = self.totalSteps + 1

            # Do the periodic work.
            self.doPeriodicWork()

        # Force the periodic work to fire at the end of an advance (except for any redistribution).
        if maxSteps != 0:
            thpt = self.redistribute
            self.redistribute = None
            self.doPeriodicWork(force=True)
            self.redistribute = thpt

        db = self.integrator.dataBase
        bcs = self.integrator.uniqueBoundaryConditions()
        numActualGhostNodes = 0
        for bc in bcs:
            numActualGhostNodes += bc.numGhostNodes
        numInternal = db.globalNumInternalNodes
        numGhost = db.globalNumGhostNodes
        numActGhost = mpi.allreduce(numActualGhostNodes, mpi.SUM)
        print(f"Total number of (internal, ghost, active ghost) nodes : ({numInternal}, {numGhost}, {numActGhost})")
        adiak_value("total_internal_nodes", numInternal)
        adiak_value("total_ghost_nodes", numGhost)
        adiak_value("total_steps", self.totalSteps)
        TimerMgr.timer_end("advance")

        # Print how much time was spent per integration cycle.
        self.stepTimer.printStatus()

        return

    #--------------------------------------------------------------------------
    # Step method, where the user can specify to take a given number of steps.
    #--------------------------------------------------------------------------
    def step(self, steps=1):
        self.advance(1e40, steps)

    #--------------------------------------------------------------------------
    # Do the periodic work.
    #--------------------------------------------------------------------------
    def doPeriodicWork(self, force=False):
        dt = self.lastDt()
        t1 = self.time()
        t0 = t1 - dt

        # Do any periodic cycle work, as determined by the number of steps.
        for method, frequency in self._periodicWork:
            if frequency is not None and (force or self.totalSteps % frequency == 0):
                method(self.totalSteps, t1, dt)

        # Do any periodic time work, as determined by the current time.
        for method, frequency in self._periodicTimeWork:
            if frequency is not None and (force or ((t0 // frequency) != (t1 // frequency))):
                method(self.totalSteps, t1, dt)

        return

    #--------------------------------------------------------------------------
    # Add a (method, frequency) tuple to the cyclic periodic work.
    # call during advance.
    #--------------------------------------------------------------------------
    def appendPeriodicWork(self, method, frequency):
        self._periodicWork.append((method, frequency))

    #--------------------------------------------------------------------------
    # Add a (method, frequency) tuple to the time based periodic work.
    #--------------------------------------------------------------------------
    def appendPeriodicTimeWork(self, method, frequency):
        self._periodicTimeWork.append((method, frequency))

    #--------------------------------------------------------------------------
    # Remove all instances of given method from the set of periodic work.
    #--------------------------------------------------------------------------
    def removePeriodicWork(self, method):
        cache = self._periodicWork[:]
        self._periodicWork = []
        for tup in cache:
            if tup[0] != method:
                self._periodicWork.append(tup)

    #--------------------------------------------------------------------------
    # Change the frequency at which the given method is called in periodic
    # work.
    #--------------------------------------------------------------------------
    def setFrequency(self, method, frequency):
        i = 0
        while (i < len(self._periodicWork) and
               self._periodicWork[i][0] != method):
            i += 1

        if i == len(self._periodicWork):
            print("Error, could not find periodic work calling ", method)
            return
        else:
            self._periodicWork[i] = (method, frequency)

    #--------------------------------------------------------------------------
    # Get the frequency at which the given method is called in periodic
    # work.
    #--------------------------------------------------------------------------
    def getFrequency(self, method):
        i = 0
        while (i < len(self._periodicWork) and
               self._periodicWork[i][0] != method):
            i += 1

        if i == len(self._periodicWork):
            raise "Error, could not find periodic work calling %s" % str(method)

        return self._periodicWork[i][1]

    #--------------------------------------------------------------------------
    # A method for printing the cycle information.
    #--------------------------------------------------------------------------
    def printCycleStatus(self, cycle, Time, dt):
        print("Cycle=%i, \tTime=%g, \tTimeStep=%g" % (cycle, Time, dt))
        return

    #--------------------------------------------------------------------------
    # Periodically update the conservation statistics.
    #--------------------------------------------------------------------------
    def updateConservation(self, cycle, Time, dt):
        self.conserve.updateHistory(cycle, Time)
        return
    
    #--------------------------------------------------------------------------
    # Periodically force garbage collection.
    #--------------------------------------------------------------------------
    def garbageCollection(self, cycle, Time, dt):
##         everything = globalFrame().f_globals
##         stuff = [x for x in everything if hasattr(everything[x], "__wards__")]
##         print "SpheralController.garbageCollection: Stuff I found with wards: ", stuff
##         for x in stuff:
##             y = everything[x]
##             for z in y.__wards__:
##                 del z
##             del y.__wards__
        gc.collect()
        return

    #--------------------------------------------------------------------------
    # Periodically drop a restart file.
    # We do a garbage collection clean up pass here since sometimes temporary
    # objects are restartable and wind up dumping data here.  Doesn't really
    # hurt anything, but it can be wasteful.
    #--------------------------------------------------------------------------
    def updateRestart(self, cycle, Time, dt):
        import gc
        while gc.collect():
            pass
        self.dropRestartFile()
        return

    #--------------------------------------------------------------------------
    # Periodically reinitialize neighbors.
    #--------------------------------------------------------------------------
    def reinitializeNeighbors(self, cycle, Time, dt):
        db = self.integrator.dataBase
        db.reinitializeNeighbors()
        return

    #--------------------------------------------------------------------------
    # Periodically redistribute the nodes between domains.
    #--------------------------------------------------------------------------
    def updateDomainDistribution(self, cycle, Time, dt):
        if self.redistribute:

            # It is *critical* that each NodeList have the same number of fields
            # registered against it on each processor, therefore we pause to
            # garbage collect here and make sure any temporaries are gone.
            import gc
            while gc.collect():
                pass

            self.redistributeTimer.start()
            self.redistribute.redistributeNodes(self.integrator.dataBase,
                                                self.integrator.uniqueBoundaryConditions())
            self.redistributeTimer.stop()
            self.redistributeTimer.printStatus()
        return

    #--------------------------------------------------------------------------
    # Find the name associated with the given object.
    #--------------------------------------------------------------------------
    def findname(thing):
        for mod in list(sys.modules.values()):
            for name, val in list(mod.__dict__.items()):
                if val is thing:
                    return name

    #--------------------------------------------------------------------------
    # Iterate over all the restartable objects and drop their state to a file.
    #--------------------------------------------------------------------------
    def dropRestartFile(self):

        if not self.restartBaseName:
            return
        # First find out if the requested directory exists.
        import os
        dire = os.path.dirname(os.path.abspath(self.restartBaseName))
        if not os.path.exists(dire):
            raise RuntimeError("Directory %s does not exist or is inaccessible." %
                               dire)

        # Now we can invoke the restart!
        import time
        start = time.time()
        fileName = self.restartBaseName + "_cycle%i" % self.totalSteps
        if self.restartFileConstructor is SidreFileIO and self.SPIOFileCountPerTimeslice is not None:
            file = self.restartFileConstructor(fileName, Create, self.SPIOFileCountPerTimeslice)
        else:
            file = self.restartFileConstructor(fileName, Create)
        RestartRegistrar.instance().dumpState(file)
        print("Wrote restart file in %0.2f seconds" % (time.time() - start))

        file.close()
        del file
        return

    #--------------------------------------------------------------------------
    # Iterate over all the restartable objects and restore their state.
    #--------------------------------------------------------------------------
    def loadRestartFile(self, restoreCycle,
                        frameDict=None):

        if not self.restartBaseName:
            return
        # Find out if the requested file exists.
        import os
        fileName = self.restartBaseName + "_cycle%i" % restoreCycle
        if self.restartFileConstructor is GzipFileIO:
            fileName += ".gz"
        if self.restartFileConstructor is SiloFileIO:
            fileName += ".silo"
        # Sidre already adds ".root" to the end of the file so we need to run the check without adding anything to the fileName
        if self.restartFileConstructor is SidreFileIO:
            if (not os.path.exists(fileName + ".root")) and (mpi.rank == 0) and (not mpi.is_fake_mpi()):
                raise RuntimeError("File %s does not exist or is inaccessible." %
                                   fileName)
        elif not os.path.exists(fileName):
            raise RuntimeError("File %s does not exist or is inaccessible." %
                               fileName)

        # Read that sucker.
        print('Reading from restart file', fileName)
        import time
        start = time.time()
        if self.restartFileConstructor is SidreFileIO and self.SPIOFileCountPerTimeslice is not None:
            file = self.restartFileConstructor(fileName, Read, self.SPIOFileCountPerTimeslice)
        else:
            file = self.restartFileConstructor(fileName, Read)
        RestartRegistrar.instance().restoreState(file)
        print("Finished: required %0.2f seconds" % (time.time() - start))

        # Reset neighboring.
        db = self.integrator.dataBase
        db.reinitializeNeighbors()
        db.updateConnectivityMap(False)

        # Do we need to force a boundary update to create ghost nodes?
        if (self.integrator.updateBoundaryFrequency > 1 and
            self.integrator.currentCycle % self.integrator.updateBoundaryFrequency != 0):
            print("Creating ghost nodes.")
            start = time.time()
            self.integrator.setGhostNodes()
            print("Finished: required %0.2f seconds" % (time.time() - start))

        file.close()
        del file
        return

    #--------------------------------------------------------------------------
    # Generate a label for our restart info.
    #--------------------------------------------------------------------------
    def label(self):
        return "Controller"

    #--------------------------------------------------------------------------
    # Store the controllers necessary data.
    #--------------------------------------------------------------------------
    def dumpState(self, file, path):
        controlPath = path + "/self"
        file.writeObject(self.totalSteps, controlPath + "/totalSteps")
        return

    #--------------------------------------------------------------------------
    # Restore the controllers state.
    #--------------------------------------------------------------------------
    def restoreState(self, file, path):
        controlPath = path + "/self"
        self.totalSteps = file.readObject(controlPath + "/totalSteps")
        return

    #--------------------------------------------------------------------------
    # Properly arrange the physics packages, including adding any needed ones.
    #--------------------------------------------------------------------------
    def organizePhysicsPackages(self, W, volumeType, facetedBoundaries):
        packages = self.integrator.physicsPackages()
        db = self.integrator.dataBase
        RKCorrections = eval("RKCorrections%s" % self.dim)
        vector_of_Physics = eval("vector_of_Physics%s" % self.dim)

        # Anyone require Voronoi cells?
        # If so we need the VoronoiCells physics package first
        voronoibcs = []
        index = -1
        for (ipack, package) in enumerate(packages):
            if package.requireVoronoiCells():
                pbcs = package.boundaryConditions
                voronoibcs += [bc for bc in pbcs if not bc in voronoibcs]
                if index == -1:
                    index = ipack

        if index >= 0:
            VC = eval("VoronoiCells" + self.dim)
            fb = eval("vector_of_FacetedVolume{}()".format(self.dim)) if facetedBoundaries is None else facetedBoundaries
            self.VoronoiCells = VC(kernelExtent = db.maxKernelExtent,
                                   facetedBoundaries = fb)
            for bc in voronoibcs:
                self.VoronoiCells.appendBoundary(bc)
            packages.insert(index, self.VoronoiCells)
            self.integrator.resetPhysicsPackages(packages)

        # Are there any packages that require reproducing kernels?
        # If so, insert the RKCorrections package prior to any RK packages
        rkorders = set()
        rkbcs = []
        needHessian = False
        rkUpdateInFinalize = False
        index = -1
        for (ipack, package) in enumerate(packages):
            ords = package.requireReproducingKernels()
            rkorders = rkorders.union(ords)
            needHessian |= package.requireReproducingKernelHessian()
            rkUpdateInFinalize |= package.updateReproducingKernelsInFinalize()
            if ords:
                pbcs = package.boundaryConditions
                rkbcs += [bc for bc in pbcs if not bc in rkbcs]
                if index == -1:
                    index = ipack
        if rkorders:
            if W is None:
                raise RuntimeError("SpheralController ERROR: the base interpolation kernel 'W' must be specified to the SpheralController when using Reproducing Kernels")
            self.RKCorrections = RKCorrections(orders = rkorders,
                                               dataBase = db,
                                               W = W,
                                               volumeType = volumeType,
                                               needHessian = needHessian,
                                               updateInFinalize = rkUpdateInFinalize)
            for bc in rkbcs:
                self.RKCorrections.appendBoundary(bc)
            if facetedBoundaries is not None:
                for b in facetedBoundaries:
                    self.RKCorrections.addFacetedBoundary(b)
            packages.insert(index, self.RKCorrections)
            self.integrator.resetPhysicsPackages(packages)
        
        # InflowOutflowBoundary can cause issues for other packages if it is not last
        if any("InflowOutflowBoundary" in type(p).__name__ for p in packages):
            packages.sort(key=lambda p: "InflowOutflowBoundary" in type(p).__name__)
            self.integrator.resetPhysicsPackages(packages)
        
        return

    #--------------------------------------------------------------------------
    # If necessary create and add a distributed boundary condition to each
    # physics package
    #
    # This method also enforces some priority among boundary conditions.
    # Current prioritization:
    #   1.  ConstantBoundaries
    #   2.  InflowOutflowBoundaries
    #   3.  ...
    #   4.  DistributedBoundary
    #--------------------------------------------------------------------------
    def insertDistributedBoundary(self, physicsPackages):

        # This is the list of boundary types that need to precede the distributed
        # boundary, since their ghost nodes need to be communicated.
        precedeDistributed = []
        if 2 in dims:
            precedeDistributed += [FacetedVolumeBoundary2d]            
        if 3 in dims:
            precedeDistributed += [FacetedVolumeBoundary3d,
                                   CylindricalBoundary,
                                   SphericalBoundary]
        for dim in dims:
            exec("""
precedeDistributed += [PeriodicBoundary%(dim)sd,
                       ConstantBoundary%(dim)sd,
                       InflowOutflowBoundary%(dim)sd]
""" % {"dim" : dim})

        # Check if this is a parallel process or not.
        if mpi.procs == 1:
            self.domainbc = None

        # If this is a parallel run, then we need to create a distributed
        # boundary condition and insert it into the list of boundaries for each physics
        # package.
        else:
            # exec("from SpheralCompiledPackages import NestedGridDistributedBoundary%s" % self.dim)
            # self.domainbc = eval("NestedGridDistributedBoundary%s.instance()" % self.dim)
            # from SpheralCompiledPackages import BoundingVolumeDistributedBoundary1d, \
            #                                    BoundingVolumeDistributedBoundary2d, \
            #                                    BoundingVolumeDistributedBoundary3d
            # self.domainbc = eval("BoundingVolumeDistributedBoundary%s.instance()" % self.dim)
            exec("from SpheralCompiledPackages import TreeDistributedBoundary%s" % self.dim)
            self.domainbc = eval("TreeDistributedBoundary%s.instance()" % self.dim)

        # Iterate over each of the physics packages, and arrange boundaries by priorities
        for package in physicsPackages:

            # Make a copy of the current set of boundary conditions for this package,
            # and assign priorities to enforce the desired order
            bcs = list(package.boundaryConditions)
            priorities = list(range(len(bcs)))
            for i, bc in enumerate(bcs):
                if isinstance(bc, eval("ConstantBoundary%s" % self.dim)):
                    priorities[i] = -2
                if isinstance(bc, eval("InflowOutflowBoundary%s" % self.dim)):
                    priorities[i] = -1
            #priorities, sortedbcs = (list(t) for t in zip(*sorted(zip(priorities, bcs))))
            sortedbcs = [x for _,x in sorted(zip(priorities, bcs), key=lambda tup: tup[0])]

            # Add the domain bc if needed
            if self.domainbc:
                sortedbcs.append(self.domainbc)

            # Reassign the package boundary conditions
            package.clearBoundaries()
            for bc in sortedbcs:
                package.appendBoundary(bc)

        # That's it.
        if not (self.domainbc is None):
            for package in physicsPackages:
                assert package.haveBoundary(self.domainbc)
        return

    #--------------------------------------------------------------------------
    # If this is a parallel run, then add a domain repartitioner to dynamically
    # balance the workload.
    #--------------------------------------------------------------------------
    def addRedistributeNodes(self, W):
        assert W.kernelExtent > 0.0
        self.redistribute = None
        self.redistributeTimer = SpheralTimer("Time for redistributing nodes.")
        if mpi.procs > 1:
            try:
                #self.redistribute = eval("ParmetisRedistributeNodes%s(W.kernelExtent)" % self.dim)
                #self.redistribute = eval("SortAndDivideRedistributeNodes%s(W.kernelExtent)" % self.dim)
                self.redistribute = eval("PeanoHilbertOrderRedistributeNodes%s(W.kernelExtent)" % self.dim)
                #self.redistribute = eval("VoronoiRedistributeNodes%s(W.kernelExtent)" % self.dim)
            except:
                print("Warning: this appears to be a parallel run, but Controller cannot construct")
                print("         dynamic redistributer.")
                pass
        return

    #---------------------------------------------------------------------------
    # Add visualization.
    #---------------------------------------------------------------------------
    def addVisualizationDumps(self, vizBaseName, vizDir, vizStep, vizTime,
                              vizFields, vizFieldLists):
        self.vizBaseName = vizBaseName
        self.vizDir = vizDir
        self.vizStep = vizStep
        self.vizTime = vizTime
        self.vizFields = vizFields
        self.vizFieldLists = vizFieldLists
        if vizStep is not None:
            self.appendPeriodicWork(self.dropViz, vizStep)
        if vizTime is not None:
            self.appendPeriodicTimeWork(self.updateViz, vizTime)
        return

    #---------------------------------------------------------------------------
    # Periodically drop viz files.
    #---------------------------------------------------------------------------
    def updateViz(self, cycle, Time, dt):
        self.dropViz(cycle, Time, dt)
        return

    #---------------------------------------------------------------------------
    # Actually does the viz dump.
    #---------------------------------------------------------------------------
    def dropViz(self,
                cycle = None,
                Time = None,
                dt = None):
        mpi.barrier()
        import time
        start = time.time()
        db = self.integrator.dataBase
        db.updateConnectivityMap(False)
        bcs = self.integrator.uniqueBoundaryConditions()
        self.integrator.setGhostNodes()
        self.vizMethod(self.integrator,
                       baseFileName = self.vizBaseName,
                       baseDirectory = self.vizDir,
                       fields = list(self.vizFields),
                       fieldLists = list(self.vizFieldLists),
                       currentTime = self.time(),
                       currentCycle = self.totalSteps,
                       dumpGhosts = self.vizGhosts,
                       dumpDerivatives = self.vizDerivs,
                       boundaries = bcs)
        print("Wrote viz file in %0.2f seconds" % (time.time() - start))
        return

    #--------------------------------------------------------------------------
    # Iteratively set the H tensors, until the desired convergence criteria
    # are met.
    #--------------------------------------------------------------------------
    def iterateIdealH(self, 
                      maxIdealHIterations = 50,
                      idealHTolerance = 1.0e-4):
        print("SpheralController: Initializing H's...")
        db = self.integrator.dataBase
        bcs = self.integrator.uniqueBoundaryConditions()

        # Find the smoothing scale method
        method = None
        for pkg in self.integrator.physicsPackages():
            if isinstance(pkg, eval(f"SmoothingScaleBase{self.dim}")):
                method = pkg
        if method is None:
            print("SpheralController::iterateIdealH no H update algorithm provided -- assuming standard SPH")
            method = eval(f"SPHSmoothingScale{self.dim}(IdealH, self.kernel)")
                
        packages = eval(f"vector_of_Physics{self.dim}()")
        if method.requireVoronoiCells():
            packages.append(self.VoronoiCells)
        packages.append(method)

        iterateIdealH = eval(f"iterateIdealH{self.dim}")
        iterateIdealH(db, packages, bcs, maxIdealHIterations, idealHTolerance, 0.0, False, False)

        return

    #---------------------------------------------------------------------------
    # Reinitialize the mass of each node such that the Voronoi mass density
    # matches the expected values.
    #---------------------------------------------------------------------------
    def voronoiInitializeMass(self):
        from generateMesh import generateLineMesh, generatePolygonalMesh, generatePolyhedralMesh
        db = self.integrator.dataBase
        nodeLists = db.fluidNodeLists
        boundaries = self.integrator.uniqueBoundaryConditions()
        method = eval("generate%sMesh" % {1 : "Line", 2 : "Polygonal", 3 : "Polyhedral"}[db.nDim])
        mesh, void = method(nodeLists,
                            boundaries = boundaries,
                            generateParallelConnectivity = False)
        for nodes in nodeLists:
            mass = nodes.mass()
            rho = nodes.massDensity()
            for i in range(nodes.numInternalNodes):
                mass[i] = rho[i]*mesh.zone(nodes, i).volume()
        return

    #---------------------------------------------------------------------------
    # Use a given pre-relaxation object (like the Voronoi hourglass control) to
    # relax a node distribution to some converged state.
    #---------------------------------------------------------------------------
    def prerelaxNodeDistribution(self,
                                 hourglass,
                                 rho,
                                 maxIterations = 100,
                                 tol = 1.0e-5):

        # What did the user pass in for rho?
        if type(rho) == type(1.0):
            rho = ConstantRho(rho)

        db = self.integrator.dataBase
        nodeLists = db.fluidNodeLists
        boundaries = self.integrator.uniqueBoundaryConditions()
        allpackages = self.integrator.physicsPackages()
        packages = eval("vector_of_Physics%id()" % db.nDim)
        packages.append(hourglass)

        # A helpful method for setting the density.
        def setRho():
            for nodes in nodeLists:
                pos = nodes.positions()
                massDensity = nodes.massDensity()
                for i in range(nodes.numInternalNodes):
                    massDensity[i] = rho(pos[i])
        setRho()

        # Iterate until we're done.
        iter = 0
        maxDisp = 2.0*tol
        while (iter < maxIterations and
               maxDisp > tol):
            iter += 1
            oldpos = db.fluidPosition
            oldpos.copyFields()
            state = eval("State%id(db, allpackages)" % db.nDim)
            derivs = eval("StateDerivatives%id(db, allpackages)" % db.nDim)
            self.integrator.initialize(state, derivs)
            state.update(derivs, 1.0, 0.0, 1.0)
            self.integrator.enforceBoundaries()
            self.integrator.applyGhostBoundaries()
            self.integrator.finalizeGhostBoundaries()
            if (self.integrator.postStateUpdate()):
                self.integrator.applyGhostBoundaries()
                self.integrator.finalizeGhostBoundaries()
            self.integrator.finalize(0.0, 0.0, db, state, derivs)

            # Check the displacements.
            maxDisp = 0.0
            newpos = db.fluidPosition
            for (oldf, newf) in zip(oldpos, newpos):
                maxDisp = max([(old - new).magnitude() for old, new in zip(oldf.internalValues(), newf.internalValues())] + [0.0])
            maxDisp = mpi.allreduce(maxDisp, mpi.MAX)
            print(" --> Iteration %i : max change = %g" % (iter, maxDisp))

        # That's about it.
        setRho()
        return

    #---------------------------------------------------------------------------
    # boolean check if we're running DEM instead of SPH
    #---------------------------------------------------------------------------
    def packagesContainDEM(self,packages):
        result = False
        for package in packages:
            if eval("isinstance(package,DEMBase{0}d)".format(self.integrator.dataBase.nDim)):
                result = True
        return result
