from PYB11Generator import *

#-------------------------------------------------------------------------------
# NodeListRegistrar
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11singleton
class NodeListRegistrar:

    # The instance attribute.  We expose this as a property of the class.
    @PYB11static
    @PYB11cppname("instance")
    @PYB11ignore
    @PYB11returnpolicy("reference")
    def getinstance(self):
        return "NodeListRegistrar<%(Dimension)s>&"
    instance = property(getinstance, doc="The static NodeListRegistrar<%(Dimension)s> instance.")


    # Attributes
    numNodeLists = PYB11property(doc="The number of NodeLists that have been created")
    numFluidNodeLists = PYB11property(doc="The number of FluidNodeLists that have been created")
    registeredNames = PYB11property(doc="The set of names for all NodeLists")
    registeredFluidNames = PYB11property(doc="The set of names for all FluidNodeLists")
    valid = PYB11property(doc="Internal consistency check")
    domainDecompositionIndependent = PYB11property("bool",
                                                   getter="domainDecompositionIndependent",
                                                   setter="domainDecompositionIndependent",
                                                   doc="Flag to force domain decomposition independent calculations -- some runtime penalty involved!")
