#ATS:test(SELF, label="Geometry inner/outer product unit tests")
# Unit tests for the inner & outer products.

import unittest
from math import *
from SpheralTestUtilities import fuzzyEqual

# What dimensions are we testing?
from spheralDimensions import spheralDimensions
dims = spheralDimensions()

from Spheral import *

# Create a global random number generator.
import random
random.seed(710)
ranrange = (-1.0, 1.0)

# Choose a default overall tolerance for comparisons
tol = 1.0e-7

#===============================================================================
# Compare two tensor'ish types to some tolerance.
#===============================================================================
def isEqual(x, y,
            tol = 1.0e-7):
    if hasattr(x, "__getitem__"):
        if len(x) != len(y):
            return False
        disc = sum([abs(xi - yi) for (xi, yi) in zip(x, y)])/len(x)
    else:
        disc = abs(x - y)
    return disc < tol

#===============================================================================
# Generate a random geometric type.
#===============================================================================
def fillRandom(Constructor):
    result = Constructor()
    ndim = Constructor.nDimensions
    nelem = Constructor.numElements
    for i in range(Constructor.numElements):
        result[i] = random.uniform(*ranrange)
    if "Sym" in Constructor.__name__:
        result = 0.5*(result + result.Transpose())
    return result

#===============================================================================
# Test class for inner product.
#===============================================================================
class TestInnerProduct(unittest.TestCase):

    #---------------------------------------------------------------------------
    # scalar . value
    #---------------------------------------------------------------------------
    def testScalarDotThing(self):
        for typestring in ("Vector%id", "Tensor%id", "SymTensor%id", "ThirdRankTensor%id"):
            for dim in dims:
                ttype = eval(typestring % dim)
                x = random.uniform(*ranrange)
                y = fillRandom(ttype)
                result = innerProduct(x, y)
                answer = ttype()
                for i in range(ttype.numElements):
                    answer[i] = x*y[i]
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch for %s: %s != %s" % (ttype.__name__, result, answer))
        return

    #---------------------------------------------------------------------------
    # value . scalar
    #---------------------------------------------------------------------------
    def testThingDotScalar(self):
        for typestring in ("Vector%id", "Tensor%id", "SymTensor%id", "ThirdRankTensor%id"):
            for dim in dims:
                ttype = eval(typestring % dim)
                x = random.uniform(*ranrange)
                y = fillRandom(ttype)
                result = innerProduct(y, x)
                answer = ttype()
                for i in range(ttype.numElements):
                    answer[i] = x*y[i]
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch for %s: %s != %s" % (ttype.__name__, result, answer))
        return

    #---------------------------------------------------------------------------
    # vector . vector
    #---------------------------------------------------------------------------
    def testVectorDotVector(self):
        for dim in dims:
            ttype = eval("Vector%id" % dim)
            x = fillRandom(ttype)
            y = fillRandom(ttype)
            result = innerProduct(x, y)
            answer = 0.0
            for i in range(dim):
                answer += x[i]*y[i]
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # tensor . vector
    #---------------------------------------------------------------------------
    def testTensorDotVector(self):
        for typestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                vtype = eval("Vector%id" % dim)
                ttype = eval(typestring % dim)
                x = fillRandom(ttype)
                y = fillRandom(vtype)
                result = innerProduct(x, y)
                answer = vtype()
                for i in range(dim):
                    for j in range(dim):
                        answer[i] += x(i,j)*y(j)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    #  vector . tensor
    #---------------------------------------------------------------------------
    def testVectorDotTensor(self):
        for typestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                vtype = eval("Vector%id" % dim)
                ttype = eval(typestring % dim)
                x = fillRandom(vtype)
                y = fillRandom(ttype)
                result = innerProduct(x, y)
                answer = vtype()
                for i in range(dim):
                    for j in range(dim):
                        answer[i] += x(j)*y(j,i)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # vector . tensor
    #---------------------------------------------------------------------------
    def testVectorDotTensor(self):
        for typestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                vtype = eval("Vector%id" % dim)
                ttype = eval(typestring % dim)
                x = fillRandom(vtype)
                y = fillRandom(ttype)
                result = innerProduct(x, y)
                answer = vtype()
                for i in range(dim):
                    for j in range(dim):
                        answer[j] += x(i)*y(i,j)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # thirdranktensor . vector
    #---------------------------------------------------------------------------
    def testThirdRankTensorDotVector(self):
        for dim in dims:
            vtype = eval("Vector%id" % dim)
            trttype = eval("ThirdRankTensor%id" % dim)
            ttype = eval("Tensor%id" % dim)
            x = fillRandom(trttype)
            y = fillRandom(vtype)
            result = innerProduct(x, y)
            answer = ttype()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        answer[dim*i + j] += x(i,j,k)*y(k)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch for %s.%s: %s != %s" % (trttype.__name__, ttype.__name__, result, answer))
        return

    #---------------------------------------------------------------------------
    # vector . thirdranktensor
    #---------------------------------------------------------------------------
    def testVectorDotThirdRankTensor(self):
        for dim in dims:
            vtype = eval("Vector%id" % dim)
            trttype = eval("ThirdRankTensor%id" % dim)
            ttype = eval("Tensor%id" % dim)
            x = fillRandom(vtype)
            y = fillRandom(trttype)
            result = innerProduct(x, y)
            answer = ttype()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        answer[dim*j + k] += x(i)*y(i,j,k)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch for %s.%s: %s != %s" % (vtype.__name__, trttype.__name__, result, answer))
        return

    #---------------------------------------------------------------------------
    # tensor . tensor
    #---------------------------------------------------------------------------
    def testTensorDotTensor(self):
        for t1typestring in ("Tensor%id", "SymTensor%id"):
            for t2typestring in ("Tensor%id", "SymTensor%id"):
                for dim in dims:
                    atype = eval("Tensor%id" % dim)
                    t1type = eval(t1typestring % dim)
                    t2type = eval(t2typestring % dim)
                    x = fillRandom(t1type)
                    y = fillRandom(t2type)
                    result = innerProduct(x, y)
                    answer = atype()
                    for i in range(dim):
                        for j in range(dim):
                            for k in range(dim):
                                answer[i*dim + j] += x(i,k)*y(k,j)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch for %s.%s: %s != %s, max disc=%s" % (t1type.__name__, t2type.__name__, result, answer, (result - answer).maxAbsElement()))
        return

    #---------------------------------------------------------------------------
    # thirdranktensor . tensor
    #---------------------------------------------------------------------------
    def testTensorDotThirdRankTensor(self):
        for ttypestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                trttype = eval("ThirdRankTensor%id" % dim)
                ttype = eval(ttypestring % dim)
                x = fillRandom(trttype)
                y = fillRandom(ttype)
                result = innerProduct(x, y)
                answer = trttype()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            for m in range(dim):
                                z = answer(i, j, k) + x(i, j, m)*y(m, k)
                                answer(i, j, k, z)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # tensor . thirdranktensor
    #---------------------------------------------------------------------------
    def testThirdRankTensorDotTensor(self):
        for ttypestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                trttype = eval("ThirdRankTensor%id" % dim)
                ttype = eval(ttypestring % dim)
                x = fillRandom(ttype)
                y = fillRandom(trttype)
                result = innerProduct(x, y)
                answer = trttype()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            for m in range(dim):
                                z = answer(i, j, k) + x(i, m)*y(m, j, k)
                                answer(i, j, k, z)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # fourthranktensor . vector
    #---------------------------------------------------------------------------
    def testFourthRankTensorDotVector(self):
        for dim in dims:
            vtype = eval("Vector%id" % dim)
            trttype = eval("ThirdRankTensor%id" % dim)
            frttype = eval("FourthRankTensor%id" % dim)
            x = fillRandom(frttype)
            y = fillRandom(vtype)
            result = innerProduct(x, y)
            answer = trttype()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            z = answer(i,j,k) + x(i,j,k,m)*y(m)
                            answer(i, j, k, z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # vector . fourthranktensor
    #---------------------------------------------------------------------------
    def testVectorDotFourthRankTensor(self):
        for dim in dims:
            vtype = eval("Vector%id" % dim)
            trttype = eval("ThirdRankTensor%id" % dim)
            frttype = eval("FourthRankTensor%id" % dim)
            x = fillRandom(vtype)
            y = fillRandom(frttype)
            result = innerProduct(x, y)
            answer = trttype()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            z = answer(i,j,k) + x(m)*y(m,i,j,k)
                            answer(i, j, k, z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # fourthranktensor . tensor
    #---------------------------------------------------------------------------
    def testFourthRankTensorDotTensor(self):
        for ttypestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                r2type = eval(ttypestring % dim)
                r4type = eval("FourthRankTensor%id" % dim)
                x = fillRandom(r4type)
                y = fillRandom(r2type)
                result = innerProduct(x, y)
                answer = r4type()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            for m in range(dim):
                                for n in range(dim):
                                    z = answer(i, j, k, n) + x(i, j, k, m)*y(m, n)
                                    answer(i, j, k, n, z)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # tensor . fourthranktensor
    #---------------------------------------------------------------------------
    def testTensorDotFourthRankTensor(self):
        for ttypestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                r2type = eval(ttypestring % dim)
                r4type = eval("FourthRankTensor%id" % dim)
                x = fillRandom(r2type)
                y = fillRandom(r4type)
                result = innerProduct(x, y)
                answer = r4type()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            for m in range(dim):
                                for n in range(dim):
                                    z = answer(i, k, m, n) + x(i, j)*y(j, k, m, n)
                                    answer(i, k, m, n, z)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # fourthranktensor . thirdranktensor
    #---------------------------------------------------------------------------
    def testFourthRankTensorDotThirdRankTensor(self):
        for dim in dims:
            r3type = eval("ThirdRankTensor%id" % dim)
            r4type = eval("FourthRankTensor%id" % dim)
            r5type = eval("FifthRankTensor%id" % dim)
            x = fillRandom(r4type)
            y = fillRandom(r3type)
            result = innerProduct(x, y)
            answer = r5type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            for n in range(dim):
                                for p in range(dim):
                                    z = answer(i, j, k, n, p) + x(i, j, k, m)*y(m, n, p)
                                    answer(i, j, k, n, p, z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    #  thirdranktensor . fourthranktensor
    #---------------------------------------------------------------------------
    def testThirdRankTensorDotFourthRankTensor(self):
        for dim in dims:
            r3type = eval("ThirdRankTensor%id" % dim)
            r4type = eval("FourthRankTensor%id" % dim)
            r5type = eval("FifthRankTensor%id" % dim)
            x = fillRandom(r3type)
            y = fillRandom(r4type)
            result = innerProduct(x, y)
            answer = r5type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            for n in range(dim):
                                for p in range(dim):
                                    z = answer(i, j, m, n, p) + x(i, j, k)*y(k, m, n, p)
                                    answer(i, j, m, n, p, z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # thirdranktensor . thirdranktensor
    #---------------------------------------------------------------------------
    def testThirdRankTensorDotThirdRankTensor(self):
        for dim in dims:
            r3type = eval("ThirdRankTensor%id" % dim)
            r4type = eval("FourthRankTensor%id" % dim)
            x = fillRandom(r3type)
            y = fillRandom(r3type)
            result = innerProduct(x, y)
            answer = r4type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            for n in range(dim):
                                z = answer(i, j, m, n) + x(i, j, k)*y(k, m, n)
                                answer(i, j, m, n, z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

#===============================================================================
# Test class for outer product.
#===============================================================================
class TestOuterProduct(unittest.TestCase):

    #---------------------------------------------------------------------------
    # scalar x value
    #---------------------------------------------------------------------------
    def testScalarOuterThing(self):
        for typestring in ("Vector%id", "Tensor%id", "SymTensor%id", "ThirdRankTensor%id"):
            for dim in dims:
                ttype = eval(typestring % dim)
                x = random.uniform(*ranrange)
                y = fillRandom(ttype)
                result = outerProduct(x, y)
                answer = ttype()
                for i in range(ttype.numElements):
                    answer[i] = x*y[i]
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch for %s: %s != %s" % (ttype.__name__, result, answer))
        return

    #---------------------------------------------------------------------------
    # value x scalar
    #---------------------------------------------------------------------------
    def testThingOuterScalar(self):
        for typestring in ("Vector%id", "Tensor%id", "SymTensor%id", "ThirdRankTensor%id"):
            for dim in dims:
                ttype = eval(typestring % dim)
                x = random.uniform(*ranrange)
                y = fillRandom(ttype)
                result = outerProduct(y, x)
                answer = ttype()
                for i in range(ttype.numElements):
                    answer[i] = x*y[i]
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch for %s: %s != %s" % (ttype.__name__, result, answer))
        return

    #---------------------------------------------------------------------------
    # vector x vector
    #---------------------------------------------------------------------------
    def testVectorOuterVector(self):
        for dim in dims:
            type = eval("Vector%id" % dim)
            ttype = eval("Tensor%id" % dim)
            x = fillRandom(type)
            y = fillRandom(type)
            result = outerProduct(x, y)
            answer = ttype()
            for i in range(dim):
                for j in range(dim):
                    answer(i, j, x[i]*y[j])
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # tensor x vector
    #---------------------------------------------------------------------------
    def testTensorOuterVector(self):
        for typestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                vtype = eval("Vector%id" % dim)
                ttype = eval(typestring % dim)
                trttype = eval("ThirdRankTensor%id" % dim)
                x = fillRandom(ttype)
                y = fillRandom(vtype)
                result = outerProduct(x, y)
                answer = trttype()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            answer(i, j, k, x(i,j)*y(k))
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    #  vector x tensor
    #---------------------------------------------------------------------------
    def testVectorOuterTensor(self):
        for typestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                vtype = eval("Vector%id" % dim)
                ttype = eval(typestring % dim)
                trttype = eval("ThirdRankTensor%id" % dim)
                x = fillRandom(vtype)
                y = fillRandom(ttype)
                result = outerProduct(x, y)
                answer = trttype()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            answer(i, j, k, x(i)*y(j,k))
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

#===============================================================================
# Test class for double inner product.
#===============================================================================
class TestDoubleInnerProduct(unittest.TestCase):

    #---------------------------------------------------------------------------
    # tensor .. tensor
    #---------------------------------------------------------------------------
    def testTensorDoubleDotTensor(self):
        for ttypestring1 in ("Tensor%id", "SymTensor%id"):
            for ttypestring2 in ("Tensor%id", "SymTensor%id"):
                for dim in dims:
                    t1type = eval(ttypestring1 % dim)
                    t2type = eval(ttypestring2 % dim)
                    x = fillRandom(t1type)
                    y = fillRandom(t2type)
                    result = innerDoubleProduct(x, y)
                    result2 = x.doubledot(y)
                    answer = 0.0
                    for i in range(dim):
                        for j in range(dim):
                            answer += x(i,j)*y(j,i)
                    self.assertTrue(abs(result - answer) < 1.0e-10, "Mismatch: %s != %s" % (result, answer))
                    self.assertTrue(abs(result2 - answer) < 1.0e-10, "Mismatch: %s != %s" % (result2, answer))
        return

    #---------------------------------------------------------------------------
    # tensor .. thirdranktensor
    #---------------------------------------------------------------------------
    def testTensorDoubleDotThirdRankTensor(self):
        for ttypestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                vtype = eval("Vector%id" % dim)
                r2type = eval(ttypestring % dim)
                r3type = eval("ThirdRankTensor%id" % dim)
                x = fillRandom(r2type)
                y = fillRandom(r3type)
                result = innerDoubleProduct(x, y)
                answer = vtype()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            answer[k] += x(i, j)*y(j, i, k)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # thirdranktensor .. tensor
    #---------------------------------------------------------------------------
    def testThirdRankTensorDoubleDotTensor(self):
        for ttypestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                vtype = eval("Vector%id" % dim)
                r2type = eval(ttypestring % dim)
                r3type = eval("ThirdRankTensor%id" % dim)
                x = fillRandom(r3type)
                y = fillRandom(r2type)
                result = innerDoubleProduct(x, y)
                answer = vtype()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            answer[i] += x(i, j, k)*y(k, j)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # thirdranktensor .. thirdranktensor
    #---------------------------------------------------------------------------
    def testThirdRankTensorDoubleDotThirdRankTensor(self):
        for dim in dims:
            r2type = eval("Tensor%id" % dim)
            r3type = eval("ThirdRankTensor%id" % dim)
            x = fillRandom(r3type)
            y = fillRandom(r3type)
            result = innerDoubleProduct(x, y)
            answer = r2type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            z = answer(i,m) + x(i,j,k)*y(k,j,m)
                            answer(i,m,z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # tensor .. fourthranktensor
    #---------------------------------------------------------------------------
    def testTensorDoubleDotFourthRankTensor(self):
        for ttypestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                ttype = eval("Tensor%id" % dim)
                r2type = eval(ttypestring % dim)
                r4type = eval("FourthRankTensor%id" % dim)
                x = fillRandom(r2type)
                y = fillRandom(r4type)
                result = innerDoubleProduct(x, y)
                answer = ttype()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            for m in range(dim):
                                z = answer(k, m) + x(i, j)*y(j, i, k, m)
                                answer(k, m, z)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # fourthranktensor .. tensor
    #---------------------------------------------------------------------------
    def testFourthRankTensorDoubleDotTensor(self):
        for ttypestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                ttype = eval("Tensor%id" % dim)
                r2type = eval(ttypestring % dim)
                r4type = eval("FourthRankTensor%id" % dim)
                x = fillRandom(r4type)
                y = fillRandom(r2type)
                result = innerDoubleProduct(x, y)
                answer = ttype()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            for m in range(dim):
                                z = answer(i, j) + x(i, j, k, m)*y(m, k)
                                answer(i, j, z)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # thirdranktensor .. fourthranktensor
    #---------------------------------------------------------------------------
    def testThirdRankTensorDoubleDotFourthRankTensor(self):
        for dim in dims:
            r3type = eval("ThirdRankTensor%id" % dim)
            r4type = eval("FourthRankTensor%id" % dim)
            x = fillRandom(r3type)
            y = fillRandom(r4type)
            result = innerDoubleProduct(x, y)
            answer = r3type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            for n in range(dim):
                                z = answer(i, m, n) + x(i, j, k)*y(k, j, m, n)
                                answer(i, m, n, z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    #  fourthranktensor .. thirdranktensor
    #---------------------------------------------------------------------------
    def testFourthRankTensorDoubleDotThirdRankTensor(self):
        for dim in dims:
            r3type = eval("ThirdRankTensor%id" % dim)
            r4type = eval("FourthRankTensor%id" % dim)
            x = fillRandom(r4type)
            y = fillRandom(r3type)
            result = innerDoubleProduct(x, y)
            answer = r3type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            for n in range(dim):
                                z = answer(i, j, n) + x(i, j, k, m)*y(m, k, n)
                                answer(i, j, n, z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    #  fourthranktensor .. fourthranktensor
    #---------------------------------------------------------------------------
    def testFourthRankTensorDoubleDotFourthRankTensor(self):
        for dim in dims:
            r4type = eval("FourthRankTensor%id" % dim)
            x = fillRandom(r4type)
            y = fillRandom(r4type)
            result = innerDoubleProduct(x, y)
            answer = r4type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            for n in range(dim):
                                for p in range(dim):
                                    z = answer(i, j, n, p) + x(i, j, k, m)*y(m, k, n, p)
                                    answer(i, j, n, p, z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # tensor .. fifthranktensor
    #---------------------------------------------------------------------------
    def testTensorDoubleDotFifthRankTensor(self):
        for ttypestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                r2type = eval(ttypestring % dim)
                r3type = eval("ThirdRankTensor%id" % dim)
                r5type = eval("FifthRankTensor%id" % dim)
                x = fillRandom(r2type)
                y = fillRandom(r5type)
                result = innerDoubleProduct(x, y)
                answer = r3type()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            for m in range(dim):
                                for n in range(dim):
                                    z = answer(k,m,n) + x(i,j)*y(j,i,k,m,n)
                                    answer(k,m,n,z)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # fifthranktensor .. tensor
    #---------------------------------------------------------------------------
    def testFifthRankTensorDoubleDotTensor(self):
        for ttypestring in ("Tensor%id", "SymTensor%id"):
            for dim in dims:
                r2type = eval(ttypestring % dim)
                r3type = eval("ThirdRankTensor%id" % dim)
                r5type = eval("FifthRankTensor%id" % dim)
                x = fillRandom(r5type)
                y = fillRandom(r2type)
                result = innerDoubleProduct(x, y)
                answer = r3type()
                for i in range(dim):
                    for j in range(dim):
                        for k in range(dim):
                            for m in range(dim):
                                for n in range(dim):
                                    z = answer(i,j,k) + x(i,j,k,m,n)*y(n,m)
                                    answer(i,j,k,z)
                self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # thirdranktensor .. fifthranktensor
    #---------------------------------------------------------------------------
    def testThirdRankTensorDoubleDotFifthRankTensor(self):
        for dim in dims:
            r3type = eval("ThirdRankTensor%id" % dim)
            r4type = eval("FourthRankTensor%id" % dim)
            r5type = eval("FifthRankTensor%id" % dim)
            x = fillRandom(r3type)
            y = fillRandom(r5type)
            result = innerDoubleProduct(x, y)
            answer = r4type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            for n in range(dim):
                                for p in range(dim):
                                    z = answer(i,m,n,p) + x(i,j,k)*y(k,j,m,n,p)
                                    answer(i,m,n,p,z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # fifthranktensor .. thirdranktensor
    #---------------------------------------------------------------------------
    def testFifthRankTensorDoubleDotThirdRankTensor(self):
        for dim in dims:
            r3type = eval("ThirdRankTensor%id" % dim)
            r4type = eval("FourthRankTensor%id" % dim)
            r5type = eval("FifthRankTensor%id" % dim)
            x = fillRandom(r5type)
            y = fillRandom(r3type)
            result = innerDoubleProduct(x, y)
            answer = r4type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            for n in range(dim):
                                for p in range(dim):
                                    z = answer(i,j,k,p) + x(i,j,k,m,n)*y(n,m,p)
                                    answer(i,j,k,p,z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # fourthranktensor .. fifthranktensor
    #---------------------------------------------------------------------------
    def FourthRankTensorDoubleDotFifthRankTensor(self):
        for dim in dims:
            r4type = eval("FourthRankTensor%id" % dim)
            r5type = eval("FifthRankTensor%id" % dim)
            x = fillRandom(r4type)
            y = fillRandom(r5type)
            result = innerDoubleProduct(x, y)
            answer = r5type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            for n in range(dim):
                                for p in range(dim):
                                    for q in range(dim):
                                        z = answer(i,j,n,p,q) + x(i,j,k,m)*y(m,k,n,p,q)
                                        answer(i,j,n,p,q,z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

    #---------------------------------------------------------------------------
    # fifthranktensor .. fourthranktensor
    #---------------------------------------------------------------------------
    def FifthRankTensorDoubleDotFourthRankTensor(self):
        for dim in dims:
            r4type = eval("FourthRankTensor%id" % dim)
            r5type = eval("FifthRankTensor%id" % dim)
            x = fillRandom(r5type)
            y = fillRandom(r4type)
            result = innerDoubleProduct(x, y)
            answer = r5type()
            for i in range(dim):
                for j in range(dim):
                    for k in range(dim):
                        for m in range(dim):
                            for n in range(dim):
                                for p in range(dim):
                                    for q in range(dim):
                                        z = answer(i,j,k,p,q) + x(i,j,k,m,n)*y(n,m,p,q)
                                        answer(i,j,k,p,q,z)
            self.assertTrue(isEqual(result, answer, tol=tol), "Mismatch: %s != %s" % (result, answer))
        return

if __name__ == "__main__":
    unittest.main()
