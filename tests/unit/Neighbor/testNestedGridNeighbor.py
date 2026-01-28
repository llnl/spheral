#ATS:test(SELF, np=1, level=100, label="NestedGridNeighbor unit tests")
from math import *
import unittest
import random

from Spheral import *
from NeighborTestBase import *


#===============================================================================
# Radom node distribution -- 1-D.
#===============================================================================
NeighborRandom1d._NeighborType = NestedGridNeighbor1d

#===============================================================================
# Radom node distribution -- 2-D.
#===============================================================================
NeighborRandom2d._NeighborType = NestedGridNeighbor2d

#===============================================================================
# Radom node distribution -- 3-D.
#===============================================================================
NeighborRandom3d._NeighborType = NestedGridNeighbor3d

#===============================================================================
# Cylindrical node distribution -- 2-D.
#===============================================================================
NeighborCylindrical2d._NeighborType = NestedGridNeighbor2d

def _create_skip_test():
    @unittest.skip("NestedGridNeighbor doesn't support ghost->ghost connectivity")
    def skip_test(self):
        pass
    return skip_test

# NestedGridNeighbor doesn't do ghost->ghost connectivity, so the overlap
# neighbor tests will choke.
for cls in [NeighborRandom1d, NeighborRandom2d, NeighborRandom3d, NeighborCylindrical2d]:
    cls.testConnectivityMapOverlapNeighbors = _create_skip_test()
    cls.testConnectivityComputeIntersection = _create_skip_test()

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
