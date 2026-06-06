import time
from Spheral3d import *
from PolyhedronFileUtilities import *

# Create a simple unit box.
verts0 = vector_of_Vector()
verts0.append(Vector(0,0,0))
verts0.append(Vector(1,0,0))
verts0.append(Vector(1,1,0))
verts0.append(Vector(0,1,0))
verts0.append(Vector(0,0,1))
verts0.append(Vector(1,0,1))
verts0.append(Vector(1,1,1))
verts0.append(Vector(0,1,1))
poly0 = Polyhedron(verts0)

# Now refine it.
t0 = time.clock()
poly1 = refinePolyhedron(poly0, 3)
t1 = time.clock()
print("Required %g seconds to refine polyhedron." % (t1 - t0))
writePolyhedronOFF(poly1, "refine_cube.off")
