#-------------------------------------------------------------------------------
# AdaptiveRefinement
#
# This class is intended as a basic adaptive refinement framework for splitting
# and combining nodes.  This class uses the fact that in the SpheralController
# you can register arbitrary methods to be called on a given frequency using
# control.appendPeriodicWork(method)
# 
# History:
#  2006-03-18 : Created by JMO.
#-------------------------------------------------------------------------------
from Spheral import *
from SpheralTestUtilities import fuzzyEqual

#-------------------------------------------------------------------------------
# Adaptive refinement is the topmost class you pass to the controller.
# The idea is you construct an AdaptiveRefinement object with your choices of
# how to select nodes for refinement and then perform the refinement operation.
#-------------------------------------------------------------------------------
class AdaptiveRefinement:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self,
                 dataBase,
                 selectionAlgorithm,
                 refinementAlgorithm,
                 spheralController):
        self.dataBase = dataBase
        self.selectionAlgorithm = selectionAlgorithm
        self.refinementAlgorithm = refinementAlgorithm
        self.spheralController = spheralController
        return

    #---------------------------------------------------------------------------
    # The method that actually does the refinement.
    #---------------------------------------------------------------------------
    def refineNodes(self, cycle, time, dt):

        # The number of new nodes we expect each refine node to be split into.
        numNewNodesPerRefineNode = self.refinementAlgorithm.numNewNodesPerRefineNode

        # Prepare the selection and refinement algorithms.
        db = self.dataBase
        self.selectionAlgorithm.prepareForSelection(db)
        self.refinementAlgorithm.prepareForRefinement(db)

        # Go over each NodeList.
        for nodeList in db.nodeLists:

            # Get the list of nodes that we want to refine on this NodeList.
            refineIDs = self.selectionAlgorithm.selectNodes(nodeList)

            # The number of nodes we expect this NodeList to have when we're done.
            finalNumNodes = (nodeList.numInternalNodes - len(refineIDs) +
                             numNewNodesPerRefineNode * len(refineIDs))

            # Allocate enough space for the new nodes we're going to create.
            numNewNodes = numNewNodesPerRefineNode * len(refineIDs)
            firstNewNode = nodeList.numInternalNodes
            nodeList.numInternalNodes += numNewNodes
            assert nodeList.numInternalNodes == firstNewNode + numNewNodes

            # Build the set of daughter particle IDs for each parent.
            daughters = [list(range(firstNewNode + i*numNewNodesPerRefineNode,
                               firstNewNode + (i + 1)*numNewNodesPerRefineNode))
                         for i in range(len(refineIDs))]
            assert len(daughters) == len(refineIDs)
            for x in daughters:
                assert len(x) == numNewNodesPerRefineNode

            # Define the new node properties.
            self.refinementAlgorithm.refineNodes(nodeList,
                                                 refineIDs,
                                                 daughters)
            # Iterativly update H.
            # Then use the mass sum to re-compute the mass density.
            if numNewNodes!=0: 
               self.spheralController.iterateIdealH()
               self.spheralController.computeMassDensity()
               nodeList.updateWeight()

            # Remove the old nodes from the NodeList.
            vecIDs = vector_of_int()
            vecIDs.extend(refineIDs)
            nodeList.deleteNodes(vecIDs)
            assert nodeList.numInternalNodes == finalNumNodes

        # That's it.
        return

#-------------------------------------------------------------------------------
# This class performs a simple 1-D refinement operation.
#-------------------------------------------------------------------------------
class SplitNodes1d:

    #---------------------------------------------------------------------------
    # Constructor.
    #---------------------------------------------------------------------------
    def __init__(self):
        self.numNewNodesPerRefineNode = 2
        return

    #---------------------------------------------------------------------------
    # Prepare for refinement by taking the spatial gradients of the mass
    # density, velocity, and specific thermal energy.
    #---------------------------------------------------------------------------
    def prepareForRefinement(self, dataBase):
        position = dataBase.fluidPosition
        weight = dataBase.fluidWeight
        mass = dataBase.fluidMass
        vel = dataBase.fluidVelocity
        rho = dataBase.fluidMassDensity
        eps = dataBase.fluidSpecificThermalEnergy
        H = dataBase.fluidHfield
        W = TableKernel1d(BSplineKernel1d(), 1000)
        self.gradrho = gradientScalar1d(rho, position, weight, mass, rho, H, W)
        self.gradvel = gradientVector1d(vel, position, weight, mass, rho, H, W)
        self.gradeps = gradientScalar1d(eps, position, weight, mass, rho, H, W)
        return

    #---------------------------------------------------------------------------
    # Split the conserved properties of the given nodes to the specified
    # daughters.
    #---------------------------------------------------------------------------
    def refineNodes(self,
                    nodeList,
                    parentList,
                    daughterList):

        # Get the fields from the NodeList.
        mass = nodeList.mass()
        position = nodeList.positions()
        velocity = nodeList.velocity()
        rho = nodeList.massDensity()
        eps = nodeList.specificThermalEnergy()
        H = nodeList.Hfield()

        gradrho = self.gradrho.fieldForNodeList(nodeList)
        gradvel = self.gradvel.fieldForNodeList(nodeList)
        gradeps = self.gradeps.fieldForNodeList(nodeList)

        # Iterate over each parent and it's daughters.
        for (parent, daughters) in zip(parentList, daughterList):

            # Check the input to make sure it's sensible.
            assert parent >= 0 and parent < nodeList.numInternalNodes
            assert len(daughters) == self.numNewNodesPerRefineNode
            for j in daughters:
                assert j >= 0 and j < nodeList.numInternalNodes

            j1 = daughters[0]
            j2 = daughters[1]

            # Compute the desired new node spacing and positions.
            xi = position[parent]
            hi = 1.0/H[parent].xx
            nPerh = nodeList.nodesPerSmoothingScale
            dx = Vector1d(0.25*hi/nPerh)
            x1 = xi - dx
            x2 = xi + dx
            position[j1] = x1
            position[j2] = x2

            # Set the intial H's.
            hj = 1.0*hi
            H[j1] = SymTensor1d(1.0/hj)
            H[j2] = SymTensor1d(1.0/hj)

            # Split the masses to try and reproduce the initial mass density profile.
            # Note we force the masses of the daughters to sum to the parent's mass.
            mi = mass[parent]
            rhoi = rho[parent]
            gradrhoi = gradrho[parent]
            rho1 = rhoi + gradrhoi.dot(x1 - xi)
            rho2 = rhoi + gradrhoi.dot(x2 - xi)
            assert rho1 > 0.0 and rho2 > 0.0
            # equal mass distribution
#            m1 = rho1/rhoi * mi
#            m2 = rho2/rhoi * mi
            msum = m1 + m2
            assert msum > 0.0
            m1 *= mi/msum
            m2 *= mi/msum
            assert fuzzyEqual(m1 + m2, mi)
            m1 = 0.5*mi
            m2 = 0.5*mi
            mass[j1] = m1
            mass[j2] = m2
            rho[j1] = rho1
            rho[j2] = rho2
            

            # Split the velocity to match the velocity profile, but preserve linear momentum.
            vi = velocity[parent]
            gradvi = gradvel[parent]
            v1 = vi + gradvi.dot(x1 - xi)
            v2 = vi + gradvi.dot(x2 - xi)
            pmomi = mi*vi.x
            pmomsum = m1*v1.x + m2*v2.x
            mul = pmomi*pmomsum/(pmomsum*pmomsum + 1e-50)
            v1 *= mul
            v2 *= mul
            velocity[j1] = Vector1d(v1)
            velocity[j2] = Vector1d(v2)

            # Split the thermal energy to match the specific thermal energy profile.
            epsi = eps[parent]
            gradepsi = gradeps[parent]
            eps1 = epsi + gradepsi.dot(x1 - xi)
            eps2 = epsi + gradepsi.dot(x2 - xi)
            
            # Choice 1: preserving the total thermal energy.
            thermEi = mi*epsi
            thermEsum = m1*eps1 + m2*eps2
#            mul = thermEi*thermEsum/(thermEsum*thermEsum + 1.0e-50)

            # Choice 2: preserving the total energy
            kinetici =0.5* mi*vi.x*vi.x
            kineticsum = 0.5*m1*v1.x*v1.x+0.5*m2*v2.x*v2.x
            mul = (kinetici+thermEi-kineticsum)/thermEsum

            # Adjust with the parameter mul
            eps1 *= mul
            eps2 *= mul
            eps[j1] = eps1
            eps[j2] = eps2

        # Update the weights.
        nodeList.updateWeight()

        return
    
