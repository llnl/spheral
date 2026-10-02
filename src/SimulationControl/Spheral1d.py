#-------------------------------------------------------------------------------
# Import Spheral objects, setting the 1-D objects as generic names.
#-------------------------------------------------------------------------------
import Spheral
for name in [x for x in Spheral.__dict__ if x.endswith("1d")]:
    exec("%s = Spheral.__dict__['%s']" % (name.replace("1d", ""), name))
from Spheral import *
FacetedVolume = Spheral.Box1d
