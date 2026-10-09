# Position redistributors.
source("testDistributeByPosition1d.py")
source("testDistributeByPosition2d.py")

# ParMETIS redistributor.
# source("testParmetisDistribute.py")

# The space filling curve redistributors.
source("testMortonOrderDistribute.py")
source("testPeanoHilbertOrderDistribute.py")

# DistributedBoundary unit tests.
source("testDistributed1d.py")
source("testDistributed2d.py")
source("testDistributed3d.py")

# Global Field reductions (allReduce) of geometric types.
source("testFieldReductions.py")
