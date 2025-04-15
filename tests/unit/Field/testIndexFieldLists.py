import unittest
from CXXTests import *

exec(compile(open("generate2Dsetup.py").read(), "generate2Dsetup.py", 'exec'))

#===============================================================================
# Main testing class.
#===============================================================================
class TestIndexFieldLists(unittest.TestCase):

    #===========================================================================
    # Setup
    #===========================================================================
    def setUp(self):

        self.ntests = 100

        # Make sure there are no ghost nodes to start with.
        for nodes in [nodes1, nodes2, nodes3]:
            nodes.numGhostNodes = 0
            nodes.neighbor().updateNodes()

        # Generate local scalar and vector FieldLists.
        self.sfield1 = ScalarField2d(nodes1, 1.0)
        self.sfield2 = ScalarField2d(nodes2, 1.0)
        self.sfield3 = ScalarField2d(nodes3, 1.0)
        for f, mul in [(self.sfield1, 1.0),
                       (self.sfield2, 10.0),
                       (self.sfield3, 100.0)]:
            for i in range(len(f)):
                f[i] = i*mul

        self.vfield1 = VectorField2d(nodes1, Vector2d(1.0, 1.0))
        self.vfield2 = VectorField2d(nodes2, Vector2d(-1.0, -1.0))
        self.vfield3 = VectorField2d(nodes3, Vector2d(4.5, 2.5))
        for f, mul in [(self.vfield1, 1.0),
                       (self.vfield2, 10.0),
                       (self.vfield3, 100.0)]:
            for i in range(len(f)):
                f[i] *= mul

        self.scalarfieldlist = ScalarFieldList2d()
        for f in [self.sfield1, self.sfield2, self.sfield3]:
            self.scalarfieldlist.appendField(f)

        self.vectorfieldlist = VectorFieldList2d()
        for f in [self.vfield1, self.vfield2, self.vfield3]:
            self.vectorfieldlist.appendField(f)

        # It's also useful to use FieldLists based on Fields that persist.
        self.scalarfieldlist2 = dataBase.fluidMassDensity
        self.vectorfieldlist2 = dataBase.fluidVelocity

        self.scalarfieldlists = [self.scalarfieldlist, self.scalarfieldlist2]
        self.vectorfieldlists = [self.vectorfieldlist, self.vectorfieldlist2]

        return

    #===========================================================================
    # Test the all node iterator (scalar).
    #===========================================================================
    def testAllNodeIteratorScalar(self):
        for fieldlist in self.scalarfieldlists:
            result = testIndexScalarFieldListByAllNodeIterators2d(dataBase,
                                                                  fieldlist)
            self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the all node iterator (vector).
    #===========================================================================
    def testAllNodeIteratorVector(self):
        for fieldlist in self.vectorfieldlists:
            result = testIndexVectorFieldListByAllNodeIterators2d(dataBase,
                                                                  fieldlist)
            self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the internal node iterator (scalar).
    #===========================================================================
    def testInternalNodeIteratorScalar(self):
        for fieldlist in self.scalarfieldlists:
            result = testIndexScalarFieldListByInternalNodeIterators2d(dataBase,
                                                                       fieldlist)
            self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the internal node iterator (vector).
    #===========================================================================
    def testInternalNodeIteratorVector(self):
        for fieldlist in self.vectorfieldlists:
            result = testIndexVectorFieldListByInternalNodeIterators2d(dataBase,
                                                                       fieldlist)
            self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the ghost node iterator (scalar).
    #===========================================================================
    def testGhostNodeIteratorScalar(self):
        bc = ReflectingBoundary2d(Plane2d(Vector2d(0.0, 0.0),
                                          Vector2d(0.0, 1.0)))
        bc.setGhostNodes(dataBase)
        for nodes in dataBase.nodeLists:
            assert nodes1.numGhostNodes > 0
        for fieldlist in self.scalarfieldlists:
            bc.applyFieldListGhostBoundary(fieldlist)
            result = testIndexScalarFieldListByGhostNodeIterators2d(dataBase,
                                                                    fieldlist)
            self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the ghost node iterator (vector).
    #===========================================================================
    def testGhostNodeIteratorVector(self):
        bc = ReflectingBoundary2d(Plane2d(Vector2d(0.0, 0.0),
                                          Vector2d(0.0, 1.0)))
        bc.setGhostNodes(dataBase)
        for nodes in dataBase.nodeLists:
            assert nodes1.numGhostNodes > 0
        for fieldlist in self.vectorfieldlists:
            bc.applyFieldListGhostBoundary(fieldlist)
            result = testIndexVectorFieldListByGhostNodeIterators2d(dataBase,
                                                                    fieldlist)
            self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the master node iterator (scalar).
    #===========================================================================
    def testMasterNodeIteratorScalar(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            for fieldlist in self.scalarfieldlists:
                result = testIndexScalarFieldListByMasterNodeIterators2d(dataBase,
                                                                         fieldlist)
                self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the master node iterator (vector).
    #===========================================================================
    def testMasterNodeIteratorVector(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            for fieldlist in self.vectorfieldlists:
                result = testIndexVectorFieldListByMasterNodeIterators2d(dataBase,
                                                                         fieldlist)
                self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the coarse node iterator (scalar).
    #===========================================================================
    def testCoarseNodeIteratorScalar(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            for fieldlist in self.scalarfieldlists:
                result = testIndexScalarFieldListByCoarseNodeIterators2d(dataBase,
                                                                         fieldlist)
                self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the coarse node iterator (vector).
    #===========================================================================
    def testCoarseNodeIteratorVector(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            for fieldlist in self.vectorfieldlists:
                result = testIndexVectorFieldListByCoarseNodeIterators2d(dataBase,
                                                                         fieldlist)
                self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the coarse node iterator cache (scalar).
    #===========================================================================
    def testCoarseNodeIteratorCacheScalar(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            for fieldlist in self.scalarfieldlists:
                result = testCacheIndexScalarFieldListByCoarseNodeIterators2d(dataBase,
                                                                              fieldlist)
                self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the coarse node iterator cache (vector).
    #===========================================================================
    def testCoarseNodeIteratorCacheVector(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            for fieldlist in self.vectorfieldlists:
                result = testCacheIndexVectorFieldListByCoarseNodeIterators2d(dataBase,
                                                                              fieldlist)
                self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the refine node iterator (scalar).
    #===========================================================================
    def testRefineNodeIteratorScalar(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            dataBase.setRefineNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            for fieldlist in self.scalarfieldlists:
                result = testIndexScalarFieldListByRefineNodeIterators2d(dataBase,
                                                                         fieldlist)
                self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the refine node iterator (vector).
    #===========================================================================
    def testRefineNodeIteratorVector(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            dataBase.setRefineNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            for fieldlist in self.vectorfieldlists:
                result = testIndexVectorFieldListByRefineNodeIterators2d(dataBase,
                                                                         fieldlist)
                self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the refine node iterator cache (scalar).
    #===========================================================================
    def testRefineNodeIteratorCacheScalar(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            for fieldlist in self.scalarfieldlists:
                result = testCacheIndexScalarFieldListByRefineNodeIterators2d(dataBase,
                                                                              fieldlist)
                self.assertTrue(result == "OK", result)
        return

    #===========================================================================
    # Test the refine node iterator cache (vector).
    #===========================================================================
    def testRefineNodeIteratorCacheVector(self):
        for i in range(self.ntests):
            nodes = g.choice([nodes1, nodes2, nodes3])
            inode = g.randint(0, nodes.numInternalNodes - 1)
            assert inode >= 0 and inode < nodes.numInternalNodes
            dataBase.setMasterNodeLists(nodes.positions()[inode],
                                        nodes.Hfield()[inode])
            for fieldlist in self.vectorfieldlists:
                result = testCacheIndexVectorFieldListByRefineNodeIterators2d(dataBase,
                                                                              fieldlist)
                self.assertTrue(result == "OK", result)
        return

if __name__ == "__main__":
    unittest.main()
