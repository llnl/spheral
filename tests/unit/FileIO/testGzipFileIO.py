#ATS:test(SELF, label="GzipFileIO unit tests")
from Spheral import *
from FileIOTestBase import *

import os
import unittest

#-------------------------------------------------------------------------------
# GzipFileIO tests.
#-------------------------------------------------------------------------------
class GzipFileIOTest(FileIOTestBase, unittest.TestCase):

    def setUp(self):
        self.n = 10 # 1000
        self.intmin = -2**24
        self.intmax = 2**24
        self.unsignedmin = 0
        self.unsignedmax = 2**32
        self.doublemin = -1e10
        self.doublemax = 1e10
        self.constructor = GzipFileIO

        # Size the NodeLists.
        nodes1d.numInternalNodes = self.n
        nodes2d.numInternalNodes = self.n
        nodes3d.numInternalNodes = self.n

        return

    def tearDown(self):
        return

    def removeFile(self, filename):
        os.remove(filename + ".gz")
        return

    # #-------------------------------------------------------------------------------
    # # Gzip can't handle things that turn into bytes (rather than strings), so we
    # # test those specially (for fail).  Once we're on Python 3 I think there's a fix
    # # for this with the gzip compress/decompress methods.
    # #-------------------------------------------------------------------------------
    # def testWriteBox(self):
    #     print("Writing FacetedVolumes currently unsupported by GzipFileIO -- skipping")
    #     return

    # def testWritePolygon(self):
    #     print("Writing FacetedVolumes currently unsupported by GzipFileIO -- skipping")
    #     return

    # def testWritePolyhedron(self):
    #     print("Writing FacetedVolumes currently unsupported by GzipFileIO -- skipping")
    #     return

    # def testVectorFacetedVolume1d(self):
    #     print("Writing FacetedVolumes currently unsupported by GzipFileIO -- skipping")
    #     return

    # def testVectorFacetedVolume2d(self):
    #     print("Writing FacetedVolumes currently unsupported by GzipFileIO -- skipping")
    #     return

    # def testVectorFacetedVolume3d(self):
    #     print("Writing FacetedVolumes currently unsupported by GzipFileIO -- skipping")
    #     return

#-------------------------------------------------------------------------------
# Run those tests.
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
