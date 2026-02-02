from Spheral import *
from SpheralTestUtilities import fuzzyEqual

import os
import random

g = random.Random(49982020438450)  # Fix random seed

WT1d = TableKernel1d(BSplineKernel1d(), 100)
eos1d = GammaLawGasMKS1d(2.0, 2.0)
nodes1d = makeFluidNodeList1d("nodes1d", eos1d)

WT2d = TableKernel2d(BSplineKernel2d(), 100)
eos2d = GammaLawGasMKS2d(2.0, 2.0)
nodes2d = makeFluidNodeList2d("nodes2d", eos2d)

WT3d = TableKernel3d(BSplineKernel3d(), 100)
eos3d = GammaLawGasMKS3d(2.0, 2.0)
nodes3d = makeFluidNodeList3d("nodes3d", eos3d)

#-------------------------------------------------------------------------------
# Generic FileIO tests.
#-------------------------------------------------------------------------------
class FileIOTestBase:

    #---------------------------------------------------------------------------
    # Generic test of writing/reading a type (assumes operator== is the
    # appropriate test).
    #---------------------------------------------------------------------------
    def boilerPlate(self, filename, path, x0, x1):
        f = self.constructor(filename, Write)
        #print("boilerPlate: Attempting to write ", x0)
        f.write(x0, path)
        f.close()
        f = self.constructor(filename, Read)
        f.read(x1, path)
        f.close()
        if x1 == x0:
            self.removeFile(filename)
        return x1 == x0
        self.assertTrue(x1 == x0,
                        "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # random Vector
    #---------------------------------------------------------------------------
    def randomVector1d(self):
        return Vector1d(random.uniform(self.doublemin, self.doublemax))
    
    def randomVector2d(self):
        return Vector2d(random.uniform(self.doublemin, self.doublemax),
                        random.uniform(self.doublemin, self.doublemax))
    
    def randomVector3d(self):
        return Vector3d(random.uniform(self.doublemin, self.doublemax),
                        random.uniform(self.doublemin, self.doublemax),
                        random.uniform(self.doublemin, self.doublemax))

    #---------------------------------------------------------------------------
    # random Tensor
    #---------------------------------------------------------------------------
    def randomTensor1d(self):
        return Tensor1d(random.uniform(self.doublemin, self.doublemax))
    
    def randomTensor2d(self):
        return Tensor2d(random.uniform(self.doublemin, self.doublemax), random.uniform(self.doublemin, self.doublemax),
                        random.uniform(self.doublemin, self.doublemax), random.uniform(self.doublemin, self.doublemax))
    
    def randomTensor3d(self):
        return Tensor3d(random.uniform(self.doublemin, self.doublemax), random.uniform(self.doublemin, self.doublemax), random.uniform(self.doublemin, self.doublemax),
                        random.uniform(self.doublemin, self.doublemax), random.uniform(self.doublemin, self.doublemax), random.uniform(self.doublemin, self.doublemax),
                        random.uniform(self.doublemin, self.doublemax), random.uniform(self.doublemin, self.doublemax), random.uniform(self.doublemin, self.doublemax))

    #---------------------------------------------------------------------------
    # random SymTensor
    #---------------------------------------------------------------------------
    def randomSymTensor1d(self):
        return SymTensor1d(random.uniform(self.doublemin, self.doublemax))
    
    def randomSymTensor2d(self):
        xx = random.uniform(self.doublemin, self.doublemax)
        xy = random.uniform(self.doublemin, self.doublemax)
        yy = random.uniform(self.doublemin, self.doublemax)
        return SymTensor2d(xx, xy, xy, yy)
    
    def randomSymTensor3d(self):
        xx = random.uniform(self.doublemin, self.doublemax)
        xy = random.uniform(self.doublemin, self.doublemax)
        xz = random.uniform(self.doublemin, self.doublemax)
        yy = random.uniform(self.doublemin, self.doublemax)
        yz = random.uniform(self.doublemin, self.doublemax)
        zz = random.uniform(self.doublemin, self.doublemax)
        return SymTensor3d(xx, xy, xz,
                           xy, yy, yz,
                           xz, yz, zz)

    #---------------------------------------------------------------------------
    # random ThirdRankTensor
    #---------------------------------------------------------------------------
    def randomThirdRankTensor(self, TRT):
        result = TRT()
        for i in range(TRT.nDimensions):
            for j in range(TRT.nDimensions):
                for k in range(TRT.nDimensions):
                    result(i, j, k, g.uniform(self.doublemin, self.doublemax))
        return result

    #---------------------------------------------------------------------------
    # random Plane
    #---------------------------------------------------------------------------
    def randomPlane1d(self):
        return Plane1d(self.randomVector1d(), self.randomVector1d().unitVector())

    def randomPlane2d(self):
        return Plane2d(self.randomVector2d(), self.randomVector2d().unitVector())

    def randomPlane3d(self):
        return Plane3d(self.randomVector3d(), self.randomVector3d().unitVector())

    #---------------------------------------------------------------------------
    # random FacetedVolume
    #---------------------------------------------------------------------------
    def randomBox(self):
        points = [self.randomVector1d(), self.randomVector1d()]
        return Box1d(vector_of_Vector1d(points))

    def randomPolygon(self):
        N = 10
        points = [self.randomVector2d() for i in range(N)]
        return Polygon(vector_of_Vector2d(points))

    def randomPolyhedron(self):
        N = 50
        points = [self.randomVector3d() for i in range(N)]
        return Polyhedron(vector_of_Vector3d(points))

    #---------------------------------------------------------------------------
    # int
    #---------------------------------------------------------------------------
    def testInt(self):
        x0 = g.randint(self.intmin, self.intmax)
        f = self.constructor("TestInt", Write)
        f.write_int(x0, "FileIOTestBase/a/b/c/d/TestInt")
        f.close()
        f = self.constructor("TestInt", Read)
        x1 = f.read_int("FileIOTestBase/a/b/c/d/TestInt")
        f.close()
        self.assertTrue(x1 == x0,
                        "%i != %i in int test" % (x1, x0))
        self.removeFile("TestInt")
        return

    #---------------------------------------------------------------------------
    # bool
    #---------------------------------------------------------------------------
    def testBool(self):
        x0 = g.choice([True, False])
        f = self.constructor("TestBool", Write)
        f.write_bool(x0, "FileIOTestBase/TestBool")
        f.close()
        f = self.constructor("TestBool", Read)
        x1 = f.read_bool("FileIOTestBase/TestBool")
        f.close()
        self.assertTrue(x1 == x0,
                        "%s != %s in bool test" % (x1, x0))
        self.removeFile("TestBool")
        return

    #---------------------------------------------------------------------------
    # float
    #---------------------------------------------------------------------------
    def testFloat(self):
        x0 = g.uniform(self.doublemin, self.doublemax)
        f = self.constructor("TestFloat", Write)
        f.write_double(x0, "FileIOTestBase/TestFloat")
        f.close()
        f = self.constructor("TestFloat", Read)
        x1 = f.read_double("FileIOTestBase/TestFloat")
        f.close()
        self.assertTrue(x1 == x0,
                        "%s != %s in float test" % (x1, x0))
        self.removeFile("TestFloat")
        return

    #---------------------------------------------------------------------------
    # string
    #---------------------------------------------------------------------------
    def testString(self):
        x0 = "abcdefg"
        f = self.constructor("TestString", Write)
        f.write_string(x0, "FileIOTestBase/TestString")
        f.close()
        f = self.constructor("TestString", Read)
        x1 = f.read_string("FileIOTestBase/TestString")
        f.close()
        self.assertTrue(x1 == x0,
                        "%s != %s in string test" % (x1, x0))
        self.removeFile("TestString")
        return

    #---------------------------------------------------------------------------
    # empty string
    #---------------------------------------------------------------------------
    def testEmptyString(self):
        x0 = ""
        f = self.constructor("TestEmptyString", Write)
        f.write_string(x0, "FileIOTestBase/TestEmptyString")
        f.close()
        f = self.constructor("TestEmptyString", Read)
        x1 = f.read_string("FileIOTestBase/TestEmptyString")
        f.close()
        self.assertTrue(x1 == x0,
                        "%s != %s in empty string test" % (x1, x0))
        self.removeFile("TestEmptyString")
        return

    #---------------------------------------------------------------------------
    # Vector1d
    #---------------------------------------------------------------------------
    def testVector1d(self):
        x0 = self.randomVector1d()
        x1 = Vector1d()
        result = self.boilerPlate("TestVector1d",
                                  "FileIOTestBase/TestVector1d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # Tensor1d
    #---------------------------------------------------------------------------
    def testTensor1d(self):
        x0 = self.randomTensor1d()
        x1 = Tensor1d()
        result = self.boilerPlate("TestTensor1d",
                                  "FileIOTestBase/TestTensor1d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # SymTensor1d
    #---------------------------------------------------------------------------
    def testSymTensor1d(self):
        x0 = self.randomSymTensor1d()
        x1 = SymTensor1d()
        result = self.boilerPlate("TestSymTensor1d",
                                  "FileIOTestBase/TestSymTensor1d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # ThirdRankTensor1d
    #---------------------------------------------------------------------------
    def testThirdRankTensor1d(self):
        x0 = self.randomThirdRankTensor(ThirdRankTensor1d)
        x1 = ThirdRankTensor1d()
        result = self.boilerPlate("TestThirdRankTensor1d",
                                  "FileIOTestBase/TestThirdRankTensor1d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # Vector2d
    #---------------------------------------------------------------------------
    def testVector2d(self):
        x0 = self.randomVector2d()
        x1 = Vector2d()
        result = self.boilerPlate("TestVector2d",
                                  "FileIOTestBase/TestVector2d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # Tensor2d
    #---------------------------------------------------------------------------
    def testTensor2d(self):
        x0 = self.randomTensor2d()
        x1 = Tensor2d()
        result = self.boilerPlate("TestTensor2d",
                                  "FileIOTestBase/TestTensor2d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # SymTensor2d
    #---------------------------------------------------------------------------
    def testSymTensor2d(self):
        x0 = self.randomSymTensor2d()
        x1 = SymTensor2d()
        result = self.boilerPlate("TestSymTensor2d",
                                  "FileIOTestBase/TestSymTensor2d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # ThirdRankTensor2d
    #---------------------------------------------------------------------------
    def testThirdRankTensor2d(self):
        x0 = self.randomThirdRankTensor(ThirdRankTensor2d)
        x1 = ThirdRankTensor2d()
        result = self.boilerPlate("TestThirdRankTensor2d",
                                  "FileIOTestBase/TestThirdRankTensor2d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # Vector3d
    #---------------------------------------------------------------------------
    def testVector3d(self):
        x0 = self.randomVector3d()
        x1 = Vector3d()
        result = self.boilerPlate("TestVector3d",
                                  "FileIOTestBase/TestVector3d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # Tensor3d
    #---------------------------------------------------------------------------
    def testTensor3d(self):
        x0 = self.randomTensor3d()
        x1 = Tensor3d()
        result = self.boilerPlate("TestTensor3d",
                                  "FileIOTestBase/TestTensor3d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # SymTensor3d
    #---------------------------------------------------------------------------
    def testSymTensor3d(self):
        x0 = self.randomSymTensor3d()
        x1 = SymTensor3d()
        result = self.boilerPlate("TestSymTensor3d",
                                  "FileIOTestBase/TestSymTensor3d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # ThirdRankTensor3d
    #---------------------------------------------------------------------------
    def testThirdRankTensor3d(self):
        x0 = self.randomThirdRankTensor(ThirdRankTensor3d)
        x1 = ThirdRankTensor3d()
        result = self.boilerPlate("TestThirdRankTensor3d",
                                  "FileIOTestBase/TestThirdRankTensor3d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # Plane1d
    #---------------------------------------------------------------------------
    def testPlane1d(self):
        x0 = self.randomPlane1d()
        x1 = Plane1d()
        result = self.boilerPlate("TestPlane1d",
                                  "FileIOTestBase/TestPlane1d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # Plane2d
    #---------------------------------------------------------------------------
    def testPlane2d(self):
        x0 = self.randomPlane2d()
        x1 = Plane2d()
        result = self.boilerPlate("TestPlane2d",
                                  "FileIOTestBase/TestPlane2d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # Plane3d
    #---------------------------------------------------------------------------
    def testPlane3d(self):
        x0 = self.randomPlane3d()
        x1 = Plane3d()
        result = self.boilerPlate("TestPlane3d",
                                  "FileIOTestBase/TestPlane3d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # FacetedVolume1d
    #---------------------------------------------------------------------------
    def testFacetedVolume1d(self):
        x0 = self.randomBox()
        x1 = Box1d()
        result = self.boilerPlate("TestFacetedVolume1d",
                                  "FileIOTestBase/TestFacetedVolume1d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # FacetedVolume2d
    #---------------------------------------------------------------------------
    def testFacetedVolume2d(self):
        x0 = self.randomPolygon()
        x1 = Polygon()
        result = self.boilerPlate("TestFacetedVolume2d",
                                  "FileIOTestBase/TestFacetedVolume2d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # FacetedVolume3d
    #---------------------------------------------------------------------------
    def testFacetedVolume3d(self):
        x0 = self.randomPolyhedron()
        x1 = Polyhedron()
        result = self.boilerPlate("TestFacetedVolume3d",
                                  "FileIOTestBase/TestFacetedVolume3d",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<int>
    #---------------------------------------------------------------------------
    def testVectorInt(self):
        for n in (0, self.n):
            x0 = vector_of_int([g.randint(self.intmin, self.intmax) for i in range(n)])
            x1 = vector_of_int()
            result = self.boilerPlate("TestVectorInt_%i" % n,
                                      "FileIOTestBase/vector_of_int",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<double>
    #---------------------------------------------------------------------------
    def testVectorDouble(self):
        for n in (0, self.n):
            x0 = vector_of_double([g.uniform(self.doublemin, self.doublemax) for i in range(n)])
            x1 = vector_of_double()
            result = self.boilerPlate("TestVectorDouble_%i" % n,
                                      "FileIOTestBase/vector_of_double",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<string>
    #---------------------------------------------------------------------------
    def testVectorString(self):
        chars = "abcdefghijklmnopqrstuvwxyz0123456789!@#$%^&*()-_=+"
        for n in (0, self.n):
            x0, x1 = vector_of_string(), vector_of_string()
            for i in range(n):
                word = ""
                wordlen = g.randint(2, 10)
                for j in range(wordlen):
                    word += g.choice(chars)
                x0.append(word)
            assert len(x0) == n
            result = self.boilerPlate("TestVectorString_%i" % n,
                                      "FileIOTestBase/vector_of_string",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<Vector1d>
    #---------------------------------------------------------------------------
    def testVectorVector1d(self):
        for n in (0, self.n):
            x0 = vector_of_Vector1d([self.randomVector1d() for i in range(n)])
            x1 = vector_of_Vector1d()
            result = self.boilerPlate("Test_vector_of_Vector1d_%i" % n,
                                      "FileIOTestBase/vector_of_Vector1d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<Vector2d>
    #---------------------------------------------------------------------------
    def testVectorVector2d(self):
        for n in (0, self.n):
            x0 = vector_of_Vector2d([self.randomVector2d() for i in range(n)])
            x1 = vector_of_Vector2d()
            result = self.boilerPlate("Test_vector_of_Vector2d_%i" % n,
                                      "FileIOTestBase/vector_of_Vector2d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<Vector3d>
    #---------------------------------------------------------------------------
    def testVectorVector3d(self):
        for n in (0, self.n):
            x0 = vector_of_Vector3d([self.randomVector3d() for i in range(n)])
            x1 = vector_of_Vector3d()
            result = self.boilerPlate("Test_vector_of_Vector3d_%i" % n,
                                      "FileIOTestBase/vector_of_Vector3d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<Tensor1d>
    #---------------------------------------------------------------------------
    def testVectorTensor1d(self):
        for n in (0, self.n):
            x0 = vector_of_Tensor1d([self.randomTensor1d() for i in range(n)])
            x1 = vector_of_Tensor1d()
            result = self.boilerPlate("Test_vector_of_Tensor1d_%i" % n,
                                      "FileIOTestBase/vector_of_Tensor1d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<Tensor2d>
    #---------------------------------------------------------------------------
    def testVectorTensor2d(self):
        for n in (0, self.n):
            x0 = vector_of_Tensor2d([self.randomTensor2d() for i in range(n)])
            x1 = vector_of_Tensor2d()
            result = self.boilerPlate("Test_vector_of_Tensor2d_%i" % n,
                                      "FileIOTestBase/vector_of_Tensor2d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<Tensor3d>
    #---------------------------------------------------------------------------
    def testVectorTensor3d(self):
        for n in (0, self.n):
            x0 = vector_of_Tensor3d([self.randomTensor3d() for i in range(n)])
            x1 = vector_of_Tensor3d()
            result = self.boilerPlate("Test_vector_of_Tensor3d_%i" % n,
                                      "FileIOTestBase/vector_of_Tensor3d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<SymTensor1d>
    #---------------------------------------------------------------------------
    def testVectorSymTensor1d(self):
        for n in (0, self.n):
            x0 = vector_of_SymTensor1d([self.randomSymTensor1d() for i in range(n)])
            x1 = vector_of_SymTensor1d()
            result = self.boilerPlate("TestVectorSymTensor1d_%i" % n,
                                      "FileIOTestBase/vector_of_SymTensor1d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<SymTensor2d>
    #---------------------------------------------------------------------------
    def testVectorSymTensor2d(self):
        for n in (0, self.n):
            x0 = vector_of_SymTensor2d([self.randomSymTensor2d() for i in range(n)])
            x1 = vector_of_SymTensor2d()
            result = self.boilerPlate("TestVectorSymTensor2d_%i" % n,
                                      "FileIOTestBase/vector_of_SymTensor2d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<SymTensor3d>
    #---------------------------------------------------------------------------
    def testVectorSymTensor3d(self):
        for n in (0, self.n):
            x0 = vector_of_SymTensor3d([self.randomSymTensor3d() for i in range(n)])
            x1 = vector_of_SymTensor3d()
            result = self.boilerPlate("TestVectorSymTensor3d_%i" % n,
                                      "FileIOTestBase/vector_of_SymTensor3d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<ThirdRankTensor1d>
    #---------------------------------------------------------------------------
    def testVectorThirdRankTensor1d(self):
        for n in (0, self.n):
            x0 = vector_of_ThirdRankTensor1d([self.randomThirdRankTensor(ThirdRankTensor1d) for i in range(n)])
            x1 = vector_of_ThirdRankTensor1d()
            result = self.boilerPlate("TestVectorThirdRankTensor1d_%i" % n,
                                      "FileIOTestBase/vector_of_ThirdRankTensor1d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<ThirdRankTensor2d>
    #---------------------------------------------------------------------------
    def testVectorThirdRankTensor2d(self):
        for n in (0, self.n):
            x0 = vector_of_ThirdRankTensor2d([self.randomThirdRankTensor(ThirdRankTensor2d) for i in range(n)])
            x1 = vector_of_ThirdRankTensor2d()
            result = self.boilerPlate("TestVectorThirdRankTensor2d_%i" % n,
                                      "FileIOTestBase/vector_of_ThirdRankTensor2d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<ThirdRankTensor3d>
    #---------------------------------------------------------------------------
    def testVectorThirdRankTensor3d(self):
        for n in (0, self.n):
            x0 = vector_of_ThirdRankTensor3d([self.randomThirdRankTensor(ThirdRankTensor3d) for i in range(n)])
            x1 = vector_of_ThirdRankTensor3d()
            result = self.boilerPlate("TestVectorThirdRankTensor3d_%i" % n,
                                      "FileIOTestBase/vector_of_ThirdRankTensor3d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<FacetedVolume1d>
    #---------------------------------------------------------------------------
    def testVectorFacetedVolume1d(self):
        for n in (0, self.n):
            x0 = vector_of_FacetedVolume1d([self.randomBox() for i in range(n)])
            x1 = vector_of_FacetedVolume1d()
            result = self.boilerPlate("TestVectorFacetedVolume1d_%i" % n,
                                      "FileIOTestBase/vector_of_FacetedVolume1d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<FacetedVolume2d>
    #---------------------------------------------------------------------------
    def testVectorFacetedVolume2d(self):
        for n in (0, self.n):
            x0 = vector_of_FacetedVolume2d([self.randomPolygon() for i in range(n)])
            x1 = vector_of_FacetedVolume2d()
            result = self.boilerPlate("TestVectorFacetedVolume2d_%i" % n,
                                      "FileIOTestBase/vector_of_FacetedVolume2d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # vector<FacetedVolume3d>
    #---------------------------------------------------------------------------
    def testVectorFacetedVolume3d(self):
        for n in (0, self.n):
            x0 = vector_of_FacetedVolume3d([self.randomPolyhedron() for i in range(n)])
            x1 = vector_of_FacetedVolume3d()
            result = self.boilerPlate("TestVectorFacetedVolume3d_%i" % n,
                                      "FileIOTestBase/vector_of_FacetedVolume3d",
                                      x0, x1)
            self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # Intfield1d
    #---------------------------------------------------------------------------
    def testIntField1d(self):
        v0 = IntField1d("int field 1d control", nodes1d)
        for i in range(self.n):
            v0[i] = g.randint(self.intmin, self.intmax)
        assert len(v0) == self.n
        f = self.constructor("TestIntField1d", Write)
        f.write(v0, "FileIOTestBase/IntField1d")
        f.close()
        f = self.constructor("TestIntField1d", Read)
        v = IntField1d("int field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/IntField1d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%i != %i @ %i of %i in IntField1d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestIntField1d")
        return

    #---------------------------------------------------------------------------
    # IntField2d
    #---------------------------------------------------------------------------
    def testIntField2d(self):
        v0 = IntField2d("int field 2d control", nodes2d)
        for i in range(self.n):
            v0[i] = g.randint(self.intmin, self.intmax)
        assert len(v0) == self.n
        f = self.constructor("TestIntField2d", Write)
        f.write(v0, "FileIOTestBase/IntField2d")
        f.close()
        f = self.constructor("TestIntField2d", Read)
        v = IntField2d("int field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/IntField2d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%i != %i @ %i of %i in IntField2d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestIntField2d")
        return

    #---------------------------------------------------------------------------
    # IntField3d
    #---------------------------------------------------------------------------
    def testIntField3d(self):
        v0 = IntField3d("int field 3d control", nodes3d)
        for i in range(self.n):
            v0[i] = g.randint(self.intmin, self.intmax)
        assert len(v0) == self.n
        f = self.constructor("TestIntField3d", Write)
        f.write(v0, "FileIOTestBase/IntField3d")
        f.close()
        f = self.constructor("TestIntField3d", Read)
        v = IntField3d("int field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/IntField3d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%i != %i @ %i of %i in IntField3d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestIntField3d")
        return

    #---------------------------------------------------------------------------
    # UnsignedField1d
    #---------------------------------------------------------------------------
    def testUnsignedField1d(self):
        v0 = UnsignedField1d("unsigned field 1d control", nodes1d)
        for i in range(self.n):
            v0[i] = g.randint(self.unsignedmin, self.unsignedmax)
        assert len(v0) == self.n
        f = self.constructor("TestUnsignedField1d", Write)
        f.write(v0, "FileIOTestBase/UnsignedField1d")
        f.close()
        f = self.constructor("TestUnsignedField1d", Read)
        v = UnsignedField1d("unsigned field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/UnsignedField1d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%i != %i @ %i of %i in UnsignedField1d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestUnsignedField1d")
        return

    #---------------------------------------------------------------------------
    # ScalarField1d
    #---------------------------------------------------------------------------
    def testScalarField1d(self):
        v0 = ScalarField1d("scalar field 1d control", nodes1d)
        for i in range(self.n):
            v0[i] = g.uniform(self.doublemin, self.doublemax)
        assert len(v0) == self.n
        f = self.constructor("TestScalarField1d", Write)
        f.write(v0, "FileIOTestBase/ScalarField1d")
        f.close()
        f = self.constructor("TestScalarField1d", Read)
        v = ScalarField1d("scalar field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/ScalarField1d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%g != %g @ %i of %i in ScalarField1d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestScalarField1d")
        return

    #---------------------------------------------------------------------------
    # ScalarField2d
    #---------------------------------------------------------------------------
    def testScalarField2d(self):
        v0 = ScalarField2d("scalar field 2d control", nodes2d)
        for i in range(self.n):
            v0[i] = g.uniform(self.doublemin, self.doublemax)
        assert len(v0) == self.n
        f = self.constructor("TestScalarField2d", Write)
        f.write(v0, "FileIOTestBase/ScalarField2d")
        f.close()
        f = self.constructor("TestScalarField2d", Read)
        v = ScalarField2d("scalar field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/ScalarField2d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%g != %g @ %i of %i in ScalarField2d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestScalarField2d")
        return

    #---------------------------------------------------------------------------
    # ScalarField3d
    #---------------------------------------------------------------------------
    def testScalarField3d(self):
        v0 = ScalarField3d("scalar field 3d control", nodes3d)
        for i in range(self.n):
            v0[i] = g.uniform(self.doublemin, self.doublemax)
        assert len(v0) == self.n
        f = self.constructor("TestScalarField3d", Write)
        f.write(v0, "FileIOTestBase/ScalarField3d")
        f.close()
        f = self.constructor("TestScalarField3d", Read)
        v = ScalarField3d("scalar field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/ScalarField3d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%g != %g @ %i of %i in ScalarField3d test" %
                            (v[i], v0[i], i, self.n))
        self.removeFile("TestScalarField3d")
        return

    #---------------------------------------------------------------------------
    # VectorField1d
    #---------------------------------------------------------------------------
    def testVectorField1d(self):
        v0 = VectorField1d("vector field 1d control", nodes1d)
        for i in range(self.n):
            v0[i] = Vector1d(g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestVectorField1d", Write)
        f.write(v0, "FileIOTestBase/VectorField1d")
        f.close()
        f = self.constructor("TestVectorField1d", Read)
        v = VectorField1d("vector field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/VectorField1d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in VectorField1d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestVectorField1d")
        return

    #---------------------------------------------------------------------------
    # VectorField2d
    #---------------------------------------------------------------------------
    def testVectorField2d(self):
        v0 = VectorField2d("vector field 2d control", nodes2d)
        for i in range(self.n):
            v0[i] = Vector2d(g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestVectorField2d", Write)
        f.write(v0, "FileIOTestBase/VectorField2d")
        f.close()
        f = self.constructor("TestVectorField2d", Read)
        v = VectorField2d("vector field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/VectorField2d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in VectorField2d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestVectorField2d")
        return

    #---------------------------------------------------------------------------
    # VectorField3d
    #---------------------------------------------------------------------------
    def testVectorField3d(self):
        v0 = VectorField3d("vector field 3d control", nodes3d)
        for i in range(self.n):
            v0[i] = Vector3d(g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestVectorField3d", Write)
        f.write(v0, "FileIOTestBase/VectorField3d")
        f.close()
        f = self.constructor("TestVectorField3d", Read)
        v = VectorField3d("vector field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/VectorField3d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in VectorField3d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestVectorField3d")
        return

    #---------------------------------------------------------------------------
    # TensorField1d
    #---------------------------------------------------------------------------
    def testTensorField1d(self):
        v0 = TensorField1d("tensor field 1d control", nodes1d)
        for i in range(self.n):
            v0[i] = Tensor1d(g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestTensorField1d", Write)
        f.write(v0, "FileIOTestBase/TensorField1d")
        f.close()
        f = self.constructor("TestTensorField1d", Read)
        v = TensorField1d("tensor field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/TensorField1d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in TensorField1d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestTensorField1d")
        return

    #---------------------------------------------------------------------------
    # TensorField2d
    #---------------------------------------------------------------------------
    def testTensorField2d(self):
        v0 = TensorField2d("tensor field 2d control", nodes2d)
        for i in range(self.n):
            v0[i] = Tensor2d(g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestTensorField2d", Write)
        f.write(v0, "FileIOTestBase/TensorField2d")
        f.close()
        f = self.constructor("TestTensorField2d", Read)
        v = TensorField2d("tensor field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/TensorField2d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in TensorField2d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestTensorField2d")
        return

    #---------------------------------------------------------------------------
    # TensorField3d
    #---------------------------------------------------------------------------
    def testTensorField3d(self):
        v0 = TensorField3d("tensor field 3d control", nodes3d)
        for i in range(self.n):
            v0[i] = Tensor3d(g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax),
                             g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestTensorField3d", Write)
        f.write(v0, "FileIOTestBase/TensorField3d")
        f.close()
        f = self.constructor("TestTensorField3d", Read)
        v = TensorField3d("tensor field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/TensorField3d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in TensorField3d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestTensorField3d")
        return

    #---------------------------------------------------------------------------
    # SymTensorField1d
    #---------------------------------------------------------------------------
    def testSymTensorField1d(self):
        v0 = SymTensorField1d("symtensor field 1d control", nodes1d)
        for i in range(self.n):
            v0[i] = SymTensor1d(g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestSymTensorField1d", Write)
        f.write(v0, "FileIOTestBase/SymTensorField1d")
        f.close()
        f = self.constructor("TestSymTensorField1d", Read)
        v = SymTensorField1d("symtensor field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/SymTensorField1d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in SymTensorField1d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestSymTensorField1d")
        return

    #---------------------------------------------------------------------------
    # SymTensorField2d
    #---------------------------------------------------------------------------
    def testSymTensorField2d(self):
        v0 = SymTensorField2d("symtensor field 2d control", nodes2d)
        for i in range(self.n):
            xx = g.uniform(self.doublemin, self.doublemax)
            xy = g.uniform(self.doublemin, self.doublemax)
            yy = g.uniform(self.doublemin, self.doublemax)
            v0[i] = SymTensor2d(xx, xy,
                                xy, yy)
        assert len(v0) == self.n
        f = self.constructor("TestSymTensorField2d", Write)
        f.write(v0, "FileIOTestBase/SymTensorField2d")
        f.close()
        f = self.constructor("TestSymTensorField2d", Read)
        v = SymTensorField2d("symtensor field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/SymTensorField2d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in SymTensorField2d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestSymTensorField2d")
        return

    #---------------------------------------------------------------------------
    # SymTensorField3d
    #---------------------------------------------------------------------------
    def testSymTensorField3d(self):
        v0 = SymTensorField3d("symtensor field 3d control", nodes3d)
        for i in range(self.n):
            xx = g.uniform(self.doublemin, self.doublemax)
            xy = g.uniform(self.doublemin, self.doublemax)
            xz = g.uniform(self.doublemin, self.doublemax)
            yy = g.uniform(self.doublemin, self.doublemax)
            yz = g.uniform(self.doublemin, self.doublemax)
            zz = g.uniform(self.doublemin, self.doublemax)
            v0[i] = SymTensor3d(xx, xy, xz,
                                xy, yy, yz,
                                xz, yz, zz)
        assert len(v0) == self.n
        f = self.constructor("TestSymTensorField3d", Write)
        f.write(v0, "FileIOTestBase/SymTensorField3d")
        f.close()
        f = self.constructor("TestSymTensorField3d", Read)
        v = SymTensorField3d("symtensor field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/SymTensorField3d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in SymTensorField3d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestSymTensorField3d")
        return

    #---------------------------------------------------------------------------
    # ThirdRankTensorField1d
    #---------------------------------------------------------------------------
    def testThirdRankTensorField1d(self):
        v0 = ThirdRankTensorField1d("third rank tensor field 1d control", nodes1d)
        for ii in range(self.n):
            v0[ii] = ThirdRankTensor1d()
            for i in range(ThirdRankTensor1d.nDimensions):
                for j in range(ThirdRankTensor1d.nDimensions):
                    for k in range(ThirdRankTensor1d.nDimensions):
                        v0[ii](i, j, k, g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestThirdRankTensorField1d", Write)
        f.write(v0, "FileIOTestBase/ThirdRankTensorField1d")
        f.close()
        f = self.constructor("TestThirdRankTensorField1d", Read)
        v = ThirdRankTensorField1d("third rank tensor field 1d test", nodes1d)
        f.read(v, "FileIOTestBase/ThirdRankTensorField1d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in ThirdRankTensorField1d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestThirdRankTensorField1d")
        return

    #---------------------------------------------------------------------------
    # ThirdRankTensorField2d
    #---------------------------------------------------------------------------
    def testThirdRankTensorField2d(self):
        v0 = ThirdRankTensorField2d("third rank tensor field 2d control", nodes2d)
        for ii in range(self.n):
            v0[ii] = ThirdRankTensor2d()
            for i in range(ThirdRankTensor2d.nDimensions):
                for j in range(ThirdRankTensor2d.nDimensions):
                    for k in range(ThirdRankTensor2d.nDimensions):
                        v0[ii](i, j, k, g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestThirdRankTensorField2d", Write)
        f.write(v0, "FileIOTestBase/ThirdRankTensorField2d")
        f.close()
        f = self.constructor("TestThirdRankTensorField2d", Read)
        v = ThirdRankTensorField2d("third rank tensor field 2d test", nodes2d)
        f.read(v, "FileIOTestBase/ThirdRankTensorField2d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in ThirdRankTensorField2d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestThirdRankTensorField2d")
        return

    #---------------------------------------------------------------------------
    # ThirdRankTensorField3d
    #---------------------------------------------------------------------------
    def testThirdRankTensorField3d(self):
        v0 = ThirdRankTensorField3d("third rank tensor field 3d control", nodes3d)
        for ii in range(self.n):
            v0[ii] = ThirdRankTensor3d()
            for i in range(ThirdRankTensor3d.nDimensions):
                for j in range(ThirdRankTensor3d.nDimensions):
                    for k in range(ThirdRankTensor3d.nDimensions):
                        v0[ii](i, j, k, g.uniform(self.doublemin, self.doublemax))
        assert len(v0) == self.n
        f = self.constructor("TestThirdRankTensorField3d", Write)
        f.write(v0, "FileIOTestBase/ThirdRankTensorField3d")
        f.close()
        f = self.constructor("TestThirdRankTensorField3d", Read)
        v = ThirdRankTensorField3d("third rank tensor field 3d test", nodes3d)
        f.read(v, "FileIOTestBase/ThirdRankTensorField3d")
        f.close()
        assert len(v) == len(v0)
        for i in range(self.n):
            self.assertTrue(v[i] == v0[i],
                            "%s != %s @ %i of %i in ThirdRankTensorField3d test" %
                            (str(v[i]), str(v0[i]), i, self.n))
        self.removeFile("TestThirdRankTensorField3d")
        return

    #---------------------------------------------------------------------------
    # IntFieldList1d
    #---------------------------------------------------------------------------
    def testIntFieldList1d(self):
        fl0 = IntFieldList1d()
        fl0.copyFields()
        fl0.appendNewField("int field 1d control", nodes1d, 0)
        for i in range(self.n):
            fl0[0][i] = g.randint(self.intmin, self.intmax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestIntFieldList1d", Write)
        f.write(fl0, "FileIOTestBase/IntFieldList1d")
        f.close()
        f = self.constructor("TestIntFieldList1d", Read)
        fl = IntFieldList1d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/IntFieldList1d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%i != %i @ %i of %i in IntFieldList1d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestIntFieldList1d")
        return

    #---------------------------------------------------------------------------
    # IntFieldList2d
    #---------------------------------------------------------------------------
    def testIntFieldList2d(self):
        fl0 = IntFieldList2d()
        fl0.copyFields()
        fl0.appendNewField("int field 2d control", nodes2d, 0)
        for i in range(self.n):
            fl0[0][i] = g.randint(self.intmin, self.intmax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestIntFieldList2d", Write)
        f.write(fl0, "FileIOTestBase/IntFieldList2d")
        f.close()
        f = self.constructor("TestIntFieldList2d", Read)
        fl = IntFieldList2d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/IntFieldList2d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%i != %i @ %i of %i in IntFieldList2d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestIntFieldList2d")
        return

    #---------------------------------------------------------------------------
    # IntFieldList3d
    #---------------------------------------------------------------------------
    def testIntFieldList3d(self):
        fl0 = IntFieldList3d()
        fl0.copyFields()
        fl0.appendNewField("int field 3d control", nodes3d, 0)
        for i in range(self.n):
            fl0[0][i] = g.randint(self.intmin, self.intmax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestIntFieldList3d", Write)
        f.write(fl0, "FileIOTestBase/IntFieldList3d")
        f.close()
        f = self.constructor("TestIntFieldList3d", Read)
        fl = IntFieldList3d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/IntFieldList3d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%i != %i @ %i of %i in IntFieldList3d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestIntFieldList3d")
        return

    #---------------------------------------------------------------------------
    # ScalarFieldList1d
    #---------------------------------------------------------------------------
    def testScalarFieldList1d(self):
        fl0 = ScalarFieldList1d()
        fl0.copyFields()
        fl0.appendNewField("scalar field 1d control", nodes1d, 0.0)
        for i in range(self.n):
            fl0[0][i] = g.uniform(self.doublemin, self.doublemax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestScalarFieldList1d", Write)
        f.write(fl0, "FileIOTestBase/ScalarFieldList1d")
        f.close()
        f = self.constructor("TestScalarFieldList1d", Read)
        fl = ScalarFieldList1d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/ScalarFieldList1d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%g != %g @ %i of %i in ScalarFieldList1d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestScalarFieldList1d")
        return

    #---------------------------------------------------------------------------
    # ScalarFieldList2d
    #---------------------------------------------------------------------------
    def testScalarFieldList2d(self):
        fl0 = ScalarFieldList2d()
        fl0.copyFields()
        fl0.appendNewField("scalar field 2d control", nodes2d, 0.0)
        for i in range(self.n):
            fl0[0][i] = g.uniform(self.doublemin, self.doublemax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestScalarFieldList2d", Write)
        f.write(fl0, "FileIOTestBase/ScalarFieldList2d")
        f.close()
        f = self.constructor("TestScalarFieldList2d", Read)
        fl = ScalarFieldList2d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/ScalarFieldList2d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%g != %g @ %i of %i in ScalarFieldList2d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestScalarFieldList2d")
        return

    #---------------------------------------------------------------------------
    # ScalarFieldList3d
    #---------------------------------------------------------------------------
    def testScalarFieldList3d(self):
        fl0 = ScalarFieldList3d()
        fl0.copyFields()
        fl0.appendNewField("scalar field 3d control", nodes3d, 0.0)
        for i in range(self.n):
            fl0[0][i] = g.uniform(self.doublemin, self.doublemax)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestScalarFieldList3d", Write)
        f.write(fl0, "FileIOTestBase/ScalarFieldList3d")
        f.close()
        f = self.constructor("TestScalarFieldList3d", Read)
        fl = ScalarFieldList3d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/ScalarFieldList3d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%g != %g @ %i of %i in ScalarFieldList3d test" %
                            (fl[0][i], fl0[0][i], i, self.n))
        self.removeFile("TestScalarFieldList3d")
        return

    #---------------------------------------------------------------------------
    # VectorFieldList1d
    #---------------------------------------------------------------------------
    def testVectorFieldList1d(self):
        fl0 = VectorFieldList1d()
        fl0.copyFields()
        fl0.appendNewField("vector field 1d control", nodes1d, Vector1d.zero)
        for i in range(self.n):
            fl0[0][i] = Vector1d(g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestVectorFieldList1d", Write)
        f.write(fl0, "FileIOTestBase/VectorFieldList1d")
        f.close()
        f = self.constructor("TestVectorFieldList1d", Read)
        fl = VectorFieldList1d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/VectorFieldList1d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl0[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in VectorFieldList1d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestVectorFieldList1d")
        return

    #---------------------------------------------------------------------------
    # VectorFieldList2d
    #---------------------------------------------------------------------------
    def testVectorFieldList2d(self):
        fl0 = VectorFieldList2d()
        fl0.copyFields()
        fl0.appendNewField("vector field 2d control", nodes2d, Vector2d.zero)
        for i in range(self.n):
            fl0[0][i] = Vector2d(g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestVectorFieldList2d", Write)
        f.write(fl0, "FileIOTestBase/VectorFieldList2d")
        f.close()
        f = self.constructor("TestVectorFieldList2d", Read)
        fl = VectorFieldList2d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/VectorFieldList2d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in VectorFieldList2d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestVectorFieldList2d")
        return

    #---------------------------------------------------------------------------
    # VectorFieldList3d
    #---------------------------------------------------------------------------
    def testVectorFieldList3d(self):
        fl0 = VectorFieldList3d()
        fl0.copyFields()
        fl0.appendNewField("vector field 3d control", nodes3d, Vector3d.zero)
        for i in range(self.n):
            fl0[0][i] = Vector3d(g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestVectorFieldList3d", Write)
        f.write(fl0, "FileIOTestBase/VectorFieldList3d")
        f.close()
        f = self.constructor("TestVectorFieldList3d", Read)
        fl = VectorFieldList3d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/VectorFieldList3d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in VectorFieldList3d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestVectorFieldList3d")
        return

    #---------------------------------------------------------------------------
    # TensorFieldList1d
    #---------------------------------------------------------------------------
    def testTensorFieldList1d(self):
        fl0 = TensorFieldList1d()
        fl0.copyFields()
        fl0.appendNewField("vector field 1d control", nodes1d, Tensor1d.zero)
        for i in range(self.n):
            fl0[0][i] = Tensor1d(g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestTensorFieldList1d", Write)
        f.write(fl0, "FileIOTestBase/TensorFieldList1d")
        f.close()
        f = self.constructor("TestTensorFieldList1d", Read)
        fl = TensorFieldList1d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/TensorFieldList1d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl0[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in TensorFieldList1d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestTensorFieldList1d")
        return

    #---------------------------------------------------------------------------
    # TensorFieldList2d
    #---------------------------------------------------------------------------
    def testTensorFieldList2d(self):
        fl0 = TensorFieldList2d()
        fl0.copyFields()
        fl0.appendNewField("vector field 2d control", nodes2d, Tensor2d.zero)
        for i in range(self.n):
            fl0[0][i] = Tensor2d(g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestTensorFieldList2d", Write)
        f.write(fl0, "FileIOTestBase/TensorFieldList2d")
        f.close()
        f = self.constructor("TestTensorFieldList2d", Read)
        fl = TensorFieldList2d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/TensorFieldList2d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in TensorFieldList2d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestTensorFieldList2d")
        return

    #---------------------------------------------------------------------------
    # TensorFieldList3d
    #---------------------------------------------------------------------------
    def testTensorFieldList3d(self):
        fl0 = TensorFieldList3d()
        fl0.copyFields()
        fl0.appendNewField("vector field 3d control", nodes3d, Tensor3d.zero)
        for i in range(self.n):
            fl0[0][i] = Tensor3d(g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax),
                                 g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestTensorFieldList3d", Write)
        f.write(fl0, "FileIOTestBase/TensorFieldList3d")
        f.close()
        f = self.constructor("TestTensorFieldList3d", Read)
        fl = TensorFieldList3d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/TensorFieldList3d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in TensorFieldList3d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestTensorFieldList3d")
        return

    #---------------------------------------------------------------------------
    # SymTensorFieldList1d
    #---------------------------------------------------------------------------
    def testSymTensorFieldList1d(self):
        fl0 = SymTensorFieldList1d()
        fl0.copyFields()
        fl0.appendNewField("vector field 1d control", nodes1d, SymTensor1d.zero)
        for i in range(self.n):
            fl0[0][i] = SymTensor1d(g.uniform(self.doublemin, self.doublemax))
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestSymTensorFieldList1d", Write)
        f.write(fl0, "FileIOTestBase/SymTensorFieldList1d")
        f.close()
        f = self.constructor("TestSymTensorFieldList1d", Read)
        fl = SymTensorFieldList1d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/SymTensorFieldList1d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl0[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in SymTensorFieldList1d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestSymTensorFieldList1d")
        return

    #---------------------------------------------------------------------------
    # SymTensorFieldList2d
    #---------------------------------------------------------------------------
    def testSymTensorFieldList2d(self):
        fl0 = SymTensorFieldList2d()
        fl0.copyFields()
        fl0.appendNewField("vector field 2d control", nodes2d, SymTensor2d.zero)
        for i in range(self.n):
            xx = g.uniform(self.doublemin, self.doublemax)
            xy = g.uniform(self.doublemin, self.doublemax)
            yy = g.uniform(self.doublemin, self.doublemax)
            fl0[0][i] = SymTensor2d(xx, xy,
                                    xy, yy)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestSymTensorFieldList2d", Write)
        f.write(fl0, "FileIOTestBase/SymTensorFieldList2d")
        f.close()
        f = self.constructor("TestSymTensorFieldList2d", Read)
        fl = SymTensorFieldList2d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/SymTensorFieldList2d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in SymTensorFieldList2d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestSymTensorFieldList2d")
        return

    #---------------------------------------------------------------------------
    # SymTensorFieldList3d
    #---------------------------------------------------------------------------
    def testSymTensorFieldList3d(self):
        fl0 = SymTensorFieldList3d()
        fl0.copyFields()
        fl0.appendNewField("vector field 3d control", nodes3d, SymTensor3d.zero)
        for i in range(self.n):
            xx = g.uniform(self.doublemin, self.doublemax)
            xy = g.uniform(self.doublemin, self.doublemax)
            xz = g.uniform(self.doublemin, self.doublemax)
            yy = g.uniform(self.doublemin, self.doublemax)
            yz = g.uniform(self.doublemin, self.doublemax)
            zz = g.uniform(self.doublemin, self.doublemax)
            fl0[0][i] = SymTensor3d(xx, xy, xz,
                                    xy, yy, yz,
                                    xz, yz, zz)
        assert len(fl0) == 1
        assert len(fl0[0]) == self.n
        f = self.constructor("TestSymTensorFieldList3d", Write)
        f.write(fl0, "FileIOTestBase/SymTensorFieldList3d")
        f.close()
        f = self.constructor("TestSymTensorFieldList3d", Read)
        fl = SymTensorFieldList3d()
        fl.copyFields()
        f.read(fl, "FileIOTestBase/SymTensorFieldList3d")
        f.close()
        assert len(fl) == len(fl0)
        assert len(fl[0]) == len(fl0[0])
        for i in range(self.n):
            self.assertTrue(fl[0][i] == fl0[0][i],
                            "%s != %s @ %i of %i in SymTensorFieldList3d test" %
                            (str(fl[0][i]), str(fl0[0][i]), i, self.n))
        self.removeFile("TestSymTensorFieldList3d")
        return

    #---------------------------------------------------------------------------
    # testWriteBox
    #---------------------------------------------------------------------------
    def testWriteBox(self):
        x0 = self.randomBox()
        x1 = Box1d()
        result = self.boilerPlate("TestBox",
                                  "FileIOTestBase/TestBox",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # testWritePolygon
    #---------------------------------------------------------------------------
    def testWritePolygon(self):
        x0 = self.randomPolygon()
        x1 = Polygon()
        result = self.boilerPlate("TestPolygon",
                                  "FileIOTestBase/TestPolygon",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # testWritePolyhedron
    #---------------------------------------------------------------------------
    def testWritePolyhedron(self):
        x0 = self.randomPolyhedron()
        x1 = Polyhedron()
        result = self.boilerPlate("TestPolyhedron", 
                                  "FileIOTestBase/TestPolyhedron",
                                  x0, x1)
        self.assertTrue(result, "FAIL: %s != %s" % (str(x1), str(x0)))

    #---------------------------------------------------------------------------
    # writeObject(int)
    #---------------------------------------------------------------------------
    def testWriteObjectInt(self):
        x0 = g.randint(self.intmin, self.intmax)
        f = self.constructor("TestInt", Write)
        f.writeObject(x0, "FileIOTestBase/TestInt")
        f.close()
        f = self.constructor("TestInt", Read)
        x1 = f.readObject("FileIOTestBase/TestInt")
        f.close()
        self.assertTrue(x1 == x0,
                        "%i != %i in int OBJECT test" % (x1, x0))
        self.removeFile("TestInt")
        return

    #---------------------------------------------------------------------------
    # writeObject(bool)
    #---------------------------------------------------------------------------
    def testWriteObjectBool(self):
        x0 = g.choice([True, False])
        f = self.constructor("TestBool", Write)
        f.writeObject(x0, "FileIOTestBase/TestBool")
        f.close()
        f = self.constructor("TestBool", Read)
        x1 = f.readObject("FileIOTestBase/TestBool")
        f.close()
        self.assertTrue(x1 == x0,
                        "%s != %s in bool OBJECT test" % (x1, x0))
        self.removeFile("TestBool")
        return

    #---------------------------------------------------------------------------
    # writeObject(string)
    #---------------------------------------------------------------------------
    def testWriteObjectString(self):
        x0 = "abcdefg"
        f = self.constructor("TestObject", Write)
        f.writeObject(x0, "FileIOTestBase/TestObject")
        f.close()
        f = self.constructor("TestObject", Read)
        x1 = f.readObject("FileIOTestBase/TestObject")
        f.close()
        self.assertTrue(x1 == x0,
                        "%s != %s in string OBJECT test" % (x1, x0))
        self.removeFile("TestObject")
        return

    #---------------------------------------------------------------------------
    # writeObject(list)
    #---------------------------------------------------------------------------
    def testWriteObjectList(self):
        x0 = [49, 492, 59392, 784761, "ackthpt"]
        f = self.constructor("TestObject", Write)
        f.writeObject(x0, "FileIOTestBase/TestObject")
        f.close()
        f = self.constructor("TestObject", Read)
        x1 = f.readObject("FileIOTestBase/TestObject")
        f.close()
        self.assertTrue(x1 == x0,
                        "%s != %s in list OBJECT test" % (x1, x0))
        self.removeFile("TestObject")
        return

    #---------------------------------------------------------------------------
    # writeObject(Vector3d)
    #---------------------------------------------------------------------------
    def testWriteObjectVector3d(self):
        x0 = Vector3d(g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax),
                      g.uniform(self.doublemin, self.doublemax))
        x1 = Vector3d()
        f = self.constructor("TestVector3d", Write)
        f.writeObject(x0, "FileIOTestBase/TestVector3d")
        f.close()
        f = self.constructor("TestVector3d", Read)
        x1 = f.readObject("FileIOTestBase/TestVector3d")
        f.close()
        self.assertTrue(x1 == x0,
                        "%s != %s in Vector3d OBJECT test" % (str(x1), str(x0)))
        self.removeFile("TestVector3d")
        return
