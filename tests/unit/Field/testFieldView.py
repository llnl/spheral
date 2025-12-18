#ATS:test(SELF, label="FieldView unit tests")

from math import *
import unittest
from Spheral2d import *
from SpheralTestUtilities import uniqueSequence

nInt = 1000
nGhost = 200
nTot = nInt + nGhost
nodes = makeVoidNodeList("test bed",
                         numInternal = nInt,
                         numGhost = nGhost)

#-------------------------------------------------------------------------------
# Set up some random generators
#-------------------------------------------------------------------------------
import random
g = random.Random()
def randScalar():
    return g.uniform(-1e10, 1e10)

class IntTraits:
    zero = 0
    def __call__(self):
        return g.randint(-10000, 100000)

class ScalarTraits:
    zero = 0.0
    def __call__(self):
        return randScalar()

class TensorTraits:
    zero = Tensor.zero
    def __call__(self):
        return Tensor(randScalar(), randScalar(),
                      randScalar(), randScalar())

class TRTTraits:
    zero = ThirdRankTensor.zero
    def __call__(self):
        result = ThirdRankTensor()
        for j in range(len(result)):
            result[j] = randScalar()
        return result

#-------------------------------------------------------------------------------
# Test class
#-------------------------------------------------------------------------------
class testFieldView(unittest.TestCase):

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def setUp(self):
        self.ntests = 100

        self.intField = IntField("integer field", nodes)
        self.scalarField = ScalarField("scalar field", nodes)
        self.tensorField = TensorField("tensor field", nodes)
        self.trtField = ThirdRankTensorField("third rank tensor field", nodes)

        self.stuff = [(self.intField,    IntTraits()),
                      (self.scalarField, ScalarTraits()),
                      (self.tensorField, TensorTraits()),
                      (self.trtField,    TRTTraits())]

        return

    #---------------------------------------------------------------------------
    # Destructor
    #---------------------------------------------------------------------------
    def tearDown(self):
        for field, generator in self.stuff:
            del field
        return

    #---------------------------------------------------------------------------
    # index (field -> span)
    #---------------------------------------------------------------------------
    def testIndexInternalField2Span(self):
        for k in range(self.ntests):
            i = g.randint(0, nTot - 1)
            for f, gen in self.stuff:
                s = f.view()
                f[i] = gen()
                assert f[i] == s[i]
        return

    #---------------------------------------------------------------------------
    # index (span -> field)
    #---------------------------------------------------------------------------
    def testIndexInternalSpan2Field(self):
        for k in range(self.ntests):
            i = g.randint(0, nTot - 1)
            for f, gen in self.stuff:
                s = f.view()
                s[i] = gen()
                assert f[i] == s[i]
        return

    #---------------------------------------------------------------------------
    # numElements
    #---------------------------------------------------------------------------
    def testNumElements(self):
        for f, gen in self.stuff:
            s = f.view()
            assert f.numElements == s.numElements == nTot
        return

    #---------------------------------------------------------------------------
    # numInternal
    #---------------------------------------------------------------------------
    def testNumInternal(self):
        for f, gen in self.stuff:
            s = f.view()
            assert f.numInternalElements == s.numInternalElements == nInt, "{} {} != {}".format(f.numInternalElements, s.numInternalElements, nInt)
        return

    #---------------------------------------------------------------------------
    # numGhosts
    #---------------------------------------------------------------------------
    def testNumGhosts(self):
        for f, gen in self.stuff:
            s = f.view()
            assert f.numGhostElements == s.numGhostElements == nGhost, "{} {} != {}".format(f.numGhostElements, s.numGhostElements, nGhost)
        return

    # #---------------------------------------------------------------------------
    # # Zero
    # #---------------------------------------------------------------------------
    # def testZero(self):
    #     for f, s, gen in self.stuff:
    #         for k in range(self.ntests):
    #             i = g.randint(0, nTot - 1)
    #             s[i] = gen()
    #         assert s != gen.zero
    #         assert f != gen.zero
    #         s.Zero()
    #         assert s == gen.zero
    #         assert f == gen.zero
    #     return

#-------------------------------------------------------------------------------
# Run the test
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
