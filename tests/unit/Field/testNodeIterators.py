import unittest
from CXXTests import *

exec(compile(open("generate2Dsetup.py").read(), "generate2Dsetup.py", 'exec'))

#===============================================================================
# Main testing class.
#===============================================================================
class TestNodeIterators(unittest.TestCase):

    def setUp(self):
        self.ntests = 100
        return

    # Test the all node iterator.
    def testAllNodeIterator(self):
        result = testGlobalAllNodeIterators2d(dataBase)
        self.assertTrue(result == "OK", result)
        return

    # Test the internal node iterator.
    def testInternalNodeIterator(self):
        result = testGlobalInternalNodeIterators2d(dataBase)
        self.assertTrue(result == "OK", result)
        return

    # Test the ghost node iterator.
    def testGhostNodeIterator(self):
        bc = ReflectingBoundary2d(Plane2d(Vector2d(0.0, 0.0),
                                          Vector2d(0.0, 1.0)))
        bc.setGhostNodes(dataBase)
        for nodes in dataBase.nodeLists:
            assert nodes1.numGhostNodes > 0
        result = testGlobalGhostNodeIterators2d(dataBase)
        self.assertTrue(result == "OK", result)
        return

    # Test the master node iterator.
    def testMasterNodeIterator(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            result = testGlobalMasterNodeIterators2d(dataBase)
            self.assertTrue(result == "OK", result)
        return

    # Test the coarse node iterator.
    def testCoarseNodeIterator(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            result = testGlobalCoarseNodeIterators2d(dataBase)
            self.assertTrue(result == "OK", result)
        return

    # Test the refine node iterator.
    def testRefineNodeIterator(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            dataBase.setRefineNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            result = testGlobalRefineNodeIterators2d(dataBase)
            self.assertTrue(result == "OK", result)
        return

if __name__ == "__main__":
    unittest.main()
