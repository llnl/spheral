# Time test to compare the GeomTensor::eigenVectors method with the Jacobi iterative
# diagonlization scheme.  GeomTensor::eigenVectors had better be winning, or we should
# just be using the Jacobi algorithm!

from Spheral import *
import random
random.seed(547957292)
from testEigen3d import randomSymTensor3d

# The number of tensors we're going to time operating on.
ntests = 10000

# Build a list of random symmetric tensors.
tensors = [randomSymTensor3d()[0] for i in range(ntests)]
assert len(tensors) == ntests

# First time the GeomTensor::eigenVectors method.
defaultTimer = SpheralTimer("GeomTensor::eigenVectors time for %i tensors." % ntests)
defaultTimer.start()
for A in tensors:
    result = A.eigenValues()
defaultTimer.stop()
defaultTimer.printStatus()

# Now time the jacobi method.
values = Vector3d()
vectors = Tensor3d()
jacobiTimer = SpheralTimer("Jacobi2() time for %i tensors." % ntests)
jacobiTimer.start()
for A in tensors:
    nrot = jacobiDiagonalize3d(A, vectors, values)
jacobiTimer.stop()
jacobiTimer.printStatus()

