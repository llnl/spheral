#ATS:test(SELF, label="Kinsol unit test")
#-------------------------------------------------------------------------------
# Setup a non-linear set of equations and solve them with Nox
#
# Solve the equations
#  x^2 + y^2 - 4 = 0
#  exp(x) + y - 1 = 0
#-------------------------------------------------------------------------------
from Spheral import *
import numpy as np
from scipy.optimize import fsolve

# Sundials currently only works for us for MPI enabled builds, so if this is a
# non-MPI version of Spheral just tap out...
import mpi
if mpi.is_fake_mpi():
    exit(0)

#-------------------------------------------------------------------------------
# Use numpy to create our reference solution we'll check against
#-------------------------------------------------------------------------------
def equations(vars):
    x, y = vars
    eq1 = x**2 + y**2 - 4.0
    eq2 = np.exp(x) + y - 1.0
    return [eq1, eq2]

# Solve using scipy
initial_guess = [1, 1]
solution = fsolve(equations, initial_guess)

print(" Numpy solution:  (x,y) = ({}, {})".format(*(solution)))
print("                 f(x,y) = ", *(equations(solution)))

#-------------------------------------------------------------------------------
# Now try using KINSOL to solve the same thing
#-------------------------------------------------------------------------------
class EquationsOperator(SolverFunction):
    def __init__(self):
        SolverFunction.__init__(self, 2)
        return

    def __call__(self, residuals, xvec):
        assert len(xvec) == len(residuals) == 2
        x, y = xvec
        eq1 = x**2 + y**2 - 4.0
        eq2 = np.exp(x) + y - 1.0
        residuals[0] = eq1
        residuals[1] = eq2
        return

solver = KINSOL()
op = EquationsOperator()
x = vector_of_double([1.0, 1.0])
solver.solve(op, x)

print("KINSOL solution:  (x,y) = ({}, {})".format(*(x)))
print("                 f(x,y) = ", *(equations(x)))

assert np.allclose(np.array(x), solution)
