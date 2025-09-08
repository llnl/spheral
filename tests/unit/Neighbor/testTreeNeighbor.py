#ATS:test(SELF, np=1, level=100, label="TreeNeighbor unit tests", nt=10)
from math import *
import unittest
import random

from Spheral import *
from NeighborTestBase import *

#===============================================================================
# Radom node distribution -- 1-D.
#===============================================================================
NeighborRandom1d._NeighborType = TreeNeighbor1d

#===============================================================================
# Radom node distribution -- 2-D.
#===============================================================================
NeighborRandom2d._NeighborType = TreeNeighbor2d

#===============================================================================
# Radom node distribution -- 3-D.
#===============================================================================
NeighborRandom3d._NeighborType = TreeNeighbor3d

#===============================================================================
# Cylindrical node distribution -- 2-D.
#===============================================================================
NeighborCylindrical2d._NeighborType = TreeNeighbor2d

#===============================================================================
# Run the tests
#===============================================================================
if __name__ == "__main__":
    unittest.main()
