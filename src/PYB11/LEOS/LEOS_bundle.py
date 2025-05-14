#-------------------------------------------------------------------------------
# LEOS_bundle
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11singleton
class LEOS_bundle:
    """LEOS_bundle

A singleton object that manages the LEOS library instance underpinning our LEOS 
equation of state objects.  Not intended for user use/intervention.
"""

    # The instance attribute.  We expose this as a property of the class.
    @PYB11static
    @PYB11cppname("instance")
    @PYB11ignore
    def getinstance(self):
        return "LEOS_bundle&"
    instance = property(getinstance, doc="The static LEOS_bundle instance.")

    # Finalize (shutdown LEOS)
    @PYB11static
    def finalize(self):
        return "void"
