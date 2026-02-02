#======================================================================
# Simple test of how quickly we can compute the eigen values/vectors
# for a set of symmetric tensors.
#======================================================================
from Spheral3d import *
from math import *
import time

import random
random.seed(941)
ranrange = 1.0e8

n = 500000
nfreq = 100
nodes = makeVoidNodeList("void", numInternal=n)
field = SymTensorField("A", nodes)
elements = []
t0 = time.clock()
for i in range(n):
    if i % nfreq == 0:
        elements = [random.uniform(-ranrange, ranrange) for i in range(6)]
    field[i] = SymTensor(elements[0], elements[1], elements[2],
                         elements[1], elements[3], elements[4],
                         elements[2], elements[4], elements[5])
vals = VectorField("B", nodes)
vecs = TensorField("C", nodes)
print("Required %g seconds to generate %i random tensors." % (time.clock() - t0, n))

t0 = time.clock()
computeEigenValues(field, vals, vecs)
print("Required %g seconds to decompose %i symmetric 3x3 tensors." % (time.clock() - t0, n))
