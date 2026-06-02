#ATS:test(SELF, label="SphericalKernel (using Biquadratic interpolation) unit tests")
from math import *
import unittest
import numpy as np

from Spheral1d import *

# Build the SphericalKernel
t0 = time.time()
W = SphericalKernel(BSplineKernel3d())
#W = SphericalBiCubicSplineKernel()
t1 = time.time()
etamax = W.etamax
print(("Required %0.4f sec to construct SphericalKernel"% (t1 - t0)))

# The r/h distance from the origin for point i (the central point we're probing around)
etavals_i = ((i + 1)*0.1 for i in range(100))

# The values of h we'll consider
hvals_i = (0.1, 0.5, 1.0, 2.0, 5.0)

#-------------------------------------------------------------------------------
# The analytic form of the quadratic bi-cubic spline from Omang et al.
#-------------------------------------------------------------------------------
def C(q):
    return q*q - 0.75*q**4 + 0.3*q**5

def D(q):
    return 2.0*(q*q - q**3) + 0.75*q**4 - 0.1*q**5

def gradC(q):
    return 2.0*q - 3.0*q**3 + 1.5*q**4

def gradD(q):
    return 4.0*q - 6.0*q**2 + 3.0*q**3 - 0.5*q**4

def sgn0(q):
    if q == 0.0:
        return 0.0
    elif q > 0.0:
        return 1.0
    else:
        return -1.0

# W()
def W3S1(rj, ri, h):
    sigj = rj/h
    sigi = ri/h
    sigdiff = abs(sigj - sigi)
    sigplus = sigj + sigi
    result = 0.0
    if sigplus <= 1.0:
        result = C(sigplus) - C(sigdiff)
    elif sigplus <= 2.0:
        if sigdiff < 1.0:
            result = -0.1 + D(sigplus) - C(sigdiff)
        else:
            result = D(sigplus) - D(sigdiff)
    else:
        if sigdiff < 1.0:
            result = 0.7 - C(sigdiff)
        elif sigdiff < 2.0:
            result = 0.8 - D(sigdiff)
    return result/(h*rj*ri)

# gradW()
def gradW3S1(rj, ri, h):
    sigj = rj/h
    sigi = ri/h
    sigdiff = abs(sigj - sigi)
    sigplus = sigj + sigi
    sgndiff = sgn0(sigi - sigj)

    if sigplus <= 1.0:
        return -W3S1(rj, ri, h)/ri + (gradC(sigplus) - gradC(sigdiff)*sgndiff)/(h*h*ri*rj)
    elif sigplus <= 2.0:
        if sigdiff < 1.0:
            return -W3S1(rj, ri, h)/ri + (gradD(sigplus) - gradC(sigdiff)*sgndiff)/(h*h*ri*rj)
        else:
            return -W3S1(rj, ri, h)/ri + (gradD(sigplus) - gradD(sigdiff)*sgndiff)/(h*h*ri*rj)
    else:
        if sigdiff < 1.0:
            return -W3S1(rj, ri, h)/ri - gradC(sigdiff)*sgndiff/(h*h*ri*rj)
        else:
            return -W3S1(rj, ri, h)/ri - gradD(sigdiff)*sgndiff/(h*h*ri*rj)

#-------------------------------------------------------------------------------
# Return a useful r_j range for a given r_i
#-------------------------------------------------------------------------------
def rprange(r, h, etastep=0.05):
    return h*np.arange(max(0.01, r/h - 2.0), r/h + 2.0, etastep)

#-------------------------------------------------------------------------------
# Test tolerances
#-------------------------------------------------------------------------------
def Wtol(etai, etaj, etamax):
    delta = abs(etai - etaj)/etamax
    if delta < 0.5:
        return 1e-6
    if delta < 0.75:
        return 1e-5
    else:
        return 0.1

def gradWtol(etai, etaj, etamax):
    delta = abs(etai - etaj)/etamax
    delta_min = min(etai, etaj)/etamax
    assert delta_min >= 0.0
    if delta_min < 0.1:
        return 2.0
    elif delta < 0.5:
        return 5e-5
    if delta < 0.75:
        return 1e-3
    else:
        return 0.1

#-------------------------------------------------------------------------------
# Error estimator
#-------------------------------------------------------------------------------
def error(val0, val1, fuzz=1.0e-8):
    return min(abs(val1 - val0), abs(val1 - val0)/max(fuzz, abs(val0)))

#-------------------------------------------------------------------------------
# Unit tests for the SphericalKernel class using an ordinary bicubic spline base
# kernel
#-------------------------------------------------------------------------------
class TestSphericalKernel(unittest.TestCase):

    def test_kernel_vals(self):
        for hi in hvals_i:
            for etai in etavals_i:
                ri = hi*etai
                for rj in rprange(ri, hi):
                    W0 = W3S1(rj, ri, hi)
                    Wij = W(Vector(rj/hi), Vector(ri/hi), 1.0/hi)
                    self.assertTrue(error(W0, Wij) < Wtol(rj/hi, ri/hi, etamax),
                                    "Kernel value outside tolerance @ (%g, %g, %g): %g != %g, %g !<= %g" % (rj, ri, hi, W0, Wij, error(W0, Wij), Wtol(rj/hi, ri/hi, etamax)))

    def test_kernel_grad_vals(self):
        for hi in hvals_i:
            for etai in etavals_i:
                ri = hi*etai
                rjfine = rprange(ri, hi, etastep=0.001)
                W3S1fine = np.array([W3S1(rj, ri, hi) for rj in rjfine])
                gradW3S1fine = np.gradient(W3S1fine, rjfine)
                for j, rj in enumerate(rprange(ri, hi, etastep=0.05)):
                    gradWij = W.grad(Vector(ri/hi), Vector(rj/hi), SymTensor(1.0/hi)).x
                    gradW0 = gradW3S1fine[50*j]
                    self.assertTrue(error(gradW0, gradWij) < gradWtol(rj/hi, ri/hi, etamax),
                                    "Kernel gradient value outside tolerance @ (%g, %g, %g): %g != %g, %g !<= %g" % (rj, ri, hi, gradW0, gradWij, error(gradW0, gradWij), gradWtol(rj/hi, ri/hi, etamax)))

    def test_simultaneous_kernel_grad_vals(self):
        self_tol = 1.0e-8
        def self_error(X0, X1, fuzz):
            return abs(X0 - X1)/max(fuzz, max(X1, X1))
        for hi in hvals_i:
            for etai in etavals_i:
                ri = hi*etai
                for rj in rprange(ri, hi, etastep=0.05):
                    Wij0 = W(Vector(rj/hi), Vector(ri/hi), 1.0/hi)
                    gradWij0 = W.grad(Vector(rj/hi), Vector(ri/hi), SymTensor(1.0/hi)).x
                    Wij1, gradWij1 = W.kernelAndGradValue(Vector(rj/hi), Vector(ri/hi), 1.0/hi)
                    self.assertTrue(self_error(Wij0, Wij1, self_tol) < self_tol,
                                    "Kernel value lookup inconsistent @ (%g, %g, %g): %g != %g, %g !<= %g" % (rj, ri, hi, Wij0, Wij1, self_error(Wij0, Wij1, self_tol), self_tol))
                    self.assertTrue(self_error(gradWij0, gradWij1.x, self_tol) < self_tol,
                                    "Kernel gradient value lookup inconsistent @ (%g, %g, %g): %g != %g, %g !<= %g" % (rj, ri, hi, gradWij0, gradWij1.x, self_error(gradWij0, gradWij1.x, self_tol), self_tol))

if __name__ == "__main__":
    unittest.main()
