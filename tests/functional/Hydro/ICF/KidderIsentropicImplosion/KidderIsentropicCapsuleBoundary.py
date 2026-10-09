#-------------------------------------------------------------------------------
# KidderIsentropicCapsuleBoundary
#
# Enforces the boundary conditions appropriate for the Kidder isentropic ICF
# capsule.
#-------------------------------------------------------------------------------
from Spheral import *
from math import *

class KidderIsentropicCapsuleBoundary1d(Boundary1d):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self, 
                 innerBoundary,  # True => inner boundary, False => outer boundary
                 integrator,     # Needed to get the time
                 answer,         # An instance of the KidderAnalytic... class
                 nodeList,       # The NodeList we're going to be applying to
                 nrGhostNodes,   # The desired radial number of ghost nodes
                 dr0,            # Initial desired spacing in r
                 ):
        Boundary1d.__init__(self)
        self.innerBoundary = innerBoundary
        self.integrator = integrator
        self.answer = answer
        self.nodes = nodeList
        self.nr = nrGhostNodes
        self.dr0 = dr0

        # Create a mapping for each type of field we're going to do boundary
        # conditions for to the method that does the job.
        self.fieldToMethodMap = {HydroFieldNames.mass : self.massBoundary,
                                 HydroFieldNames.velocity : self.velocityBoundary,
                                 HydroFieldNames.massDensity : self.rhoBoundary,
                                 HydroFieldNames.specificThermalEnergy : self.noopBoundary,
                                 HydroFieldNames.pressure : self.pressureBoundary,
                                 HydroFieldNames.soundSpeed : self.noopBoundary,
                                 HydroFieldNames.omegaGradh : self.omegaBoundary,
                                 HydroFieldNames.specificThermalEnergy + "0" : self.noopBoundary,
                                 HydroFieldNames.A_CRKSPH : self.noopBoundary,
                                 HydroFieldNames.B_CRKSPH : self.noopBoundary,
                                 HydroFieldNames.C_CRKSPH : self.noopBoundary,
                                 HydroFieldNames.gradA_CRKSPH : self.noopBoundary,
                                 HydroFieldNames.gradB_CRKSPH : self.noopBoundary,
                                 }

        return

    #---------------------------------------------------------------------------
    # Create our ghost nodes.
    #---------------------------------------------------------------------------
    def setGhostNodes(self, nodeList):
        self.addNodeList(self.nodes)
        self.nodeIDs = list(range(self.nodes.numNodes,
                             self.nodes.numNodes + self.nr))
        nodeList.numGhostNodes = self.nodes.numGhostNodes + self.nr
        boundNodes = self.accessBoundaryNodes(nodeList)
        for i in self.nodeIDs:
            boundNodes.ghostNodes().append(i)
            boundNodes.controlNodes().append(0)

        self.updateGhostNodes(nodeList)
        return

    #---------------------------------------------------------------------------
    # Update the minimal ghost node info (positions and H's).
    #---------------------------------------------------------------------------
    def updateGhostNodes(self, nodeList):

        # Fields we're going to update.
        pos = nodeList.positions()
        H = nodeList.Hfield()
        t = self.integrator.currentTime

        # Initial values for r and H.
        if self.innerBoundary:
            dr = -self.dr0 * self.answer.hfrac(t)
            ri = self.answer.rInner(t) - 0.5*dr
        else:
            dr = self.dr0 * self.answer.hfrac(t)
            ri = self.answer.rOuter(t) - 0.5*dr
        self.dr = abs(dr)
        Hi = SymTensor1d(1.0/(self.dr*nodeList.nodesPerSmoothingScale))

        # Step out from the Boundary incrementally, updating each ghost node.
        for i in self.nodeIDs:
            ri += dr
            pos[i].x = ri
            H[i] = Hi

        mass = nodeList.mass()
        self.massBoundary(mass)
        print("mass: ", [mass[i] for i in self.nodeIDs])

        nodeList.neighbor().updateNodes()
        return

    #---------------------------------------------------------------------------
    # Update ghost field values.
    # This is specialized for each of the fields that gets passed in for a 
    # hydro run.
    #---------------------------------------------------------------------------
    def applyGhostBoundary(self, field):

        if field.name not in self.fieldToMethodMap:
            print("KidderBoundary WARNING: unable to apply boundary condition to %s" % field.name)
            return

        print("Apply GHOST to ", field.name)
        self.fieldToMethodMap[field.name](field)
        return

    #---------------------------------------------------------------------------
    # Violation nodes -- not defined in this case.
    #---------------------------------------------------------------------------
    def setViolationNodes(self, nodeList):
        return

    def updateViolationNodes(self, nodeList):
        return

    def enforceBoundary(self, field):
        return

    #---------------------------------------------------------------------------
    # No op boundary -- do nothing.
    #---------------------------------------------------------------------------
    def noopBoundary(self, field):
        return

    #---------------------------------------------------------------------------
    # Mass.
    #---------------------------------------------------------------------------
    def massBoundary(self, field):
        assert field.name == HydroFieldNames.mass
        t = self.integrator.currentTime
        pos = self.nodes.positions()
        print("Applying mass boundary")
        for i in self.nodeIDs:
            rhoi = self.answer.rho(t, pos[i].x)
            field[i] = rhoi * self.dr
            assert field[i] > 0.0
            print(" --> ", i, rhoi, field[i])
        return

    #---------------------------------------------------------------------------
    # Velocity.
    #---------------------------------------------------------------------------
    def velocityBoundary(self, field):
        assert field.name == HydroFieldNames.velocity
        t = self.integrator.currentTime
        pos = self.nodes.positions()
        for i in self.nodeIDs:
            field[i].x = self.answer.vr(t, pos[i].x)
        return

    #---------------------------------------------------------------------------
    # Mass density.
    #---------------------------------------------------------------------------
    def rhoBoundary(self, field):
        assert field.name == HydroFieldNames.massDensity
        t = self.integrator.currentTime
        pos = self.nodes.positions()
        for i in self.nodeIDs:
            field[i] = self.answer.rho(t, pos[i].x)
        return

    #---------------------------------------------------------------------------
    # Pressure.
    #---------------------------------------------------------------------------
    def pressureBoundary(self, field):
        assert field.name == HydroFieldNames.pressure
        t = self.integrator.currentTime
        pos = self.nodes.positions()
        for i in self.nodeIDs:
            field[i] = self.answer.P(t, pos[i].x)
        return

    #---------------------------------------------------------------------------
    # Omega.
    #---------------------------------------------------------------------------
    def omegaBoundary(self, field):
        assert field.name == HydroFieldNames.omegaGradh
        for i in self.nodeIDs:
            field[i] = 1.0
        return

#-------------------------------------------------------------------------------
# This class directly sets the state of a set of internal nodes, rather than
# create ghost nodes.  This is necessary for the inner boundary, where ghost
# nodes would often need a negative mass density.
#-------------------------------------------------------------------------------
class KidderIsentropicCapsuleEnforcementBoundary1d(Physics1d):

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self, 
                 integrator,     # Needed to get the time
                 answer,         # An instance of the KidderAnalytic... class
                 nodeList,       # The NodeList we're going to be applying to
                 nodeIDs,        # The set of internal nodes we're going to control
                 interiorNodeIDs,# The interior nodes not controlled by a boundary condition
                 hinitial,       # Initial smoothing scale that we're going to impose
                 hmin,           # minimum smoothing scale
                 hmax,           # maximum smoothing scale
                 dr0,            # Intial node spacing
                 ):
        Physics1d.__init__(self)
        self.integrator = integrator
        self.answer = answer
        self.nodes = nodeList
        self.nodeIDs = nodeIDs
        self.interiorNodeIDs = interiorNodeIDs
        self.hinitial = hinitial
        self.hmin = hmin
        self.hmax = hmax
        self.dr0 = dr0
        self.gamma = self.answer.gamma
        self.gamma1 = self.answer.gamma1

        # Store some of the initial values for nodes we're controlling.
        pos = self.nodes.positions()
        eps = self.nodes.specificThermalEnergy()
        self.initialRadii = [pos[i].x for i in self.nodeIDs]
        self.initialEps = [eps[i] for i in self.nodeIDs]

        # Store the total mass in interior nodes.
        mass = self.nodes.mass()
        self.interiorMass = mpi.allreduce(sum([mass[i] for i in interiorNodeIDs]), mpi.SUM)

        return

    #---------------------------------------------------------------------------
    # Override the state values for the physics variables, imposing the analytic
    # solution.
    #---------------------------------------------------------------------------
    def overrideState(self, state, t):

        # Extract the state we're going to set.
        rho = state.scalarFields(HydroFieldNames.massDensity)
        vol = state.scalarFields(HydroFieldNames.volume)
        pos = state.vectorFields(HydroFieldNames.position)
        eps = state.scalarFields(HydroFieldNames.specificThermalEnergy)
        vel = state.vectorFields(HydroFieldNames.velocity)
        H = state.symTensorFields(HydroFieldNames.H)
        P = state.scalarFields(HydroFieldNames.pressure)
        cs = state.scalarFields(HydroFieldNames.soundSpeed)
        mass = state.scalarFields(HydroFieldNames.mass)
        
        hfrac = self.answer.hfrac(t)
        Hi = SymTensor1d(1.0/max(self.hmin, min(self.hmax, self.hinitial*hfrac)))

        # Now set these variables for the nodes we're controlling.
        for k in range(len(self.nodeIDs)):
            i = self.nodeIDs[k]
            ri = self.initialRadii[k]*hfrac
            voli = self.dr0*hfrac
            rhoi = mass[0][i]/voli   # self.answer.rho(t, ri)
            epsi = self.initialEps[k]/(hfrac*hfrac)
            Pi = self.gamma1*rhoi*epsi
            #Pi = self.answer.P(t, ri)
            #epsi = Pi/(self.gamma1*rhoi)
            if vol.numFields > 0:
                vol[0][i] = voli
            rho[0][i] = rhoi
            pos[0][i].x = ri
            eps[0][i] = epsi
            vel[0][i].x = self.answer.vr(t, ri)
            H[0][i] = Hi
            P[0][i] = Pi
            cs[0][i] = sqrt(self.gamma * self.gamma1 * epsi)
        
        n = self.nodes.neighbor()
        n.updateNodes()        

        return

    #---------------------------------------------------------------------------
    # Override the derivative values for the physics variables, imposing the
    # analytic solution.
    #---------------------------------------------------------------------------
    def overrideDerivatives(self, derivs, t):

        # Extract the state we're going to set.
        DxDt = derivs.vectorFields("delta " + HydroFieldNames.position)
        DvDt = derivs.vectorFields(HydroFieldNames.hydroAcceleration)
        DvDx = derivs.tensorFields(HydroFieldNames.velocityGradient)
        internalDvDx = derivs.tensorFields(HydroFieldNames.internalVelocityGradient)
        DHDt = derivs.symTensorFields("delta " + HydroFieldNames.H)
        DrhoDt = derivs.scalarFields("delta " + HydroFieldNames.massDensity)
        DepsDt = derivs.scalarFields("delta " + HydroFieldNames.specificThermalEnergy)
        #rhoSum = derivs.scalarFields("new " + HydroFieldNames.massDensity)
        Hideal = derivs.symTensorFields("new " + HydroFieldNames.H)
        XSPHDeltaV = derivs.vectorFields(HydroFieldNames.XSPHDeltaV)
        mass = self.nodes.mass()

        hfrac = self.answer.hfrac(t)
        hfracdot = self.answer.hfracDot(t)
        Hi = SymTensor1d(1.0/max(self.hmin, min(self.hmax, self.hinitial*hfrac)))

        # Now set these variables for the nodes we're controlling.
        for k in range(len(self.nodeIDs)):
            i = self.nodeIDs[k]
            ri = self.initialRadii[k]*hfrac
            rhoi = mass[i]/(self.dr0*hfrac)   # self.answer.rho(t, ri)
            epsi = self.initialEps[k]/(hfrac*hfrac)
            Pi = self.gamma1*rhoi*epsi
            #Pi = self.answer.P(t, ri)
            DvDxi = Tensor1d(self.answer.DvrDr(t, ri))
            rhoDoti = self.answer.rhoDot(t, ri)
            Pdoti = self.answer.Pdot(t, ri)

            DxDt[0][i].x = self.answer.vr(t, ri)
            DvDt[0][i].x = self.answer.vrDot(t, ri)
            DvDx[0][i] = DvDxi
            internalDvDx[0][i] = DvDxi
            DHDt[0][i] = Hi*rhoDoti/rhoi
            DrhoDt[0][i] = rhoDoti
            DepsDt[0][i] = -2.0*epsi*hfracdot/hfrac
            #DepsDt[0][i] = (Pdoti - Pi*rhoDoti/rhoi)/(self.gamma1*rhoi)
            #rhoSum[0][i] = rhoi
            Hideal[0][i] = Hi
            XSPHDeltaV.Zero()
        
        return

    #---------------------------------------------------------------------------
    # Physics::evaluateDerivatives
    #---------------------------------------------------------------------------
    def evaluateDerivatives(self, currentTime, dt, dataBase, state, derivs):
        self.overrideDerivatives(derivs, currentTime)
        return

    #---------------------------------------------------------------------------
    # Physics::dt
    #---------------------------------------------------------------------------
    def dt(self, dataBase, state, derivs, currentTime):
        return (1.0e50, "No vote.")

    #---------------------------------------------------------------------------
    # Physics::registerState
    #---------------------------------------------------------------------------
    def registerState(self, dataBase, state):
        return

    #---------------------------------------------------------------------------
    # Physics::registerDerivatives
    #---------------------------------------------------------------------------
    def registerDerivatives(self, dataBase, derivs):
        return

    #---------------------------------------------------------------------------
    # Physics::initialize
    #---------------------------------------------------------------------------
    def initialize(self, currentTime, dt, dataBase, state, derivs):
        #self.overrideState(state, currentTime)
        #self.overrideDerivatives(derivs, currentTime)
        return

    #---------------------------------------------------------------------------
    # Physics::finalize
    #---------------------------------------------------------------------------
    def finalize(self, currentTime, dt, dataBase, state, derivs):
        #self.overrideState(state, currentTime)
        return

    #---------------------------------------------------------------------------
    # Physics::finalizeDerivatives
    #---------------------------------------------------------------------------
    def finalizeDerivatives(self, currentTime, dt, dataBase, state, derivs):
        self.overrideDerivatives(derivs, currentTime)
        return

    #---------------------------------------------------------------------------
    # Physics::postStateUpdate
    #---------------------------------------------------------------------------
    def postStateUpdate(self, currentTime, dt, dataBase, state, derivs):

##         # Take a snapshot of the internal energies.
##         Key = pair_NodeList1d_string
##         mKey = Key(self.nodes, HydroFieldNames.mass)
##         epsKey = Key(self.nodes, HydroFieldNames.specificThermalEnergy)
##         mass = state.scalarField(mKey)
##         eps = state.scalarField(epsKey)
##         self.preEps = [eps[i] for i in self.nodeIDs]

        self.overrideState(state, currentTime)

##         # Redistribute the energy change we imposed to the interior nodes.
##         deltaE = sum([mass[i] * (eps[i] - self.preEps[k]) for i, k in zip(self.nodeIDs, range(len(self.nodeIDs)))])
##         print "Energy changed by ", deltaE
##         deltaEps = -deltaE/self.interiorMass
##         for i in self.interiorNodeIDs:
##             eps[i] += deltaEps

    #---------------------------------------------------------------------------
    # Physics::enforceBoundaries
    #---------------------------------------------------------------------------
    def enforceBoundaries(self, state, derivs):

        t = self.integrator.currentTime
        self.overrideState(state, t)

        return
