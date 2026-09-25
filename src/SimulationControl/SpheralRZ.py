#-------------------------------------------------------------------------------
# Import Spheral objects, setting the 2D objects as generic names.
# This version is specialized for the 2D RZ formulation.
#-------------------------------------------------------------------------------
import Spheral
for name in [x for x in Spheral.__dict__ if x.endswith("2d")]:
    exec("%s = Spheral.__dict__['%s']" % (name.replace("2d", ""), name))
for name in [x for x in Spheral.__dict__ if x != "RZ" and x.endswith("RZ")]:
    exec("%s = Spheral.__dict__['%s']" % (name.replace("RZ", ""), name))
from Spheral import *
GeometryRegistrar.coords(CoordinateType.RZ)
