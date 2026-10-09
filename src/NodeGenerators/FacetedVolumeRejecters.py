import mpi
from RejecterBase import RejecterBase
from SpheralCompiledPackages import Vector2d, Vector3d

#-------------------------------------------------------------------------------
# Reject in a polygonal surface
#-------------------------------------------------------------------------------
class PolygonalSurfaceRejecter(RejecterBase):

    def __init__(self, surface,
                 interior = True):
        RejecterBase.__init__(self)
        self.surface = surface
        self.interior = interior # Reject interior to surface?
        return

    # New style rejecter
    def accept(self, x, y):
        return self.interior ^ self.surface.contains(Vector2d(x,y))

#-------------------------------------------------------------------------------
# Reject in a polyhedral surface
#-------------------------------------------------------------------------------
class PolyhedralSurfaceRejecter(RejecterBase):

    def __init__(self, surface,
                 interior = True):
        RejecterBase.__init__(self)
        self.surface = surface
        self.interior = interior # Reject interior to surface?
        return

    # New style rejecter
    def accept(self, x, y, z):
        return self.interior ^ self.surface.contains(Vector3d(x,y,z))

