from Spheral1d import *
from math import *
import numpy as np

#-------------------------------------------------------------------------------
# The analytic answer for the acoustic wave problem.
#-------------------------------------------------------------------------------
class AcousticWaveSolution:

    # Constructor.
    def __init__(self,
                 eos,
                 nPoints=101,
                 cs=1.0,
                 rho0=1.0,
                 x0=0.0,
                 x1=1.0,
                 A=1.e-6,
                 k=1.0,
                 h0=1.0):
        self.eos=eos
        self.nPoints = nPoints
        self.cs = cs
        self.rho0 = rho0
        self.x0 = x0
        self.x1 = x1
        self.A = A
        self.k = k
        self.h0 = h0
        return

    # Compute and return the solution on the given positions.
    def solution(self, t, xvals=None):
        if xvals is None:
            xvals = np.linspace(self.x0, self.x1, self.nPoints)

        omegat = self.k*self.cs*t
        length = self.x1 - self.x0
        v = self.cs * self.A * np.sin(self.k*(xvals - self.x0)/length - omegat)
        u = np.full(len(xvals), 1.0)
        rho = self.rho0 * (1.0 + self.A * np.sin(self.k*(xvals - self.x0)/length - omegat))
        h = self.h0 * self.rho0/rho
        P = np.array([self.eos.pressure(rhoi,ui) for (rhoi,ui) in zip(rho,u)])
        return xvals, v, u, rho, P, h
