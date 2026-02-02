#-------------------------------------------------------------------------------
# Compute the analytic solution to the piston driven shock compaction of
# porosity.  Based on the 1D compaction wave problem (5.1) in 
# Jutzi, M., Benz, W., & Michel, P. (2008). Numerical simulations of impacts
# involving porous bodies.I. Implementing sub-resolution porosity in a 3D SPH
# hydrocode. Icarus, 198(1), 242–255.
#
# We are using the solution method described more fully in the section on
# Hugoniots and Isentropes from
# Herrmann, W. (1969). Constitutive Equation for the Dynamic Compaction of 
# Ductile Porous Materials. Journal of Applied Physics, 40(6), 2490–2499.
#-------------------------------------------------------------------------------
import mpi
import numpy as np
from numpy import linalg as LA
import scipy.integrate, scipy.optimize, scipy.stats
from computeHugoniotWithPorosity import *

#-------------------------------------------------------------------------------
# Compute the solution profiles for the test
#-------------------------------------------------------------------------------
class PlanarCompactionSolution:

    def __init__(self,
                 eos,
                 vpiston,
                 eps0,
                 alpha0,
                 alphat,
                 Pe,
                 Pt,
                 Ps,
                 n1,
                 n2,
                 cS0 = None,
                 c0 = None,
                 h0 = None,
                 nPoints = 101,
                 pistonFrame = False):
        self.eos = eos
        self.vpiston = abs(vpiston)
        self.eps0 = eps0
        self.alpha0 = alpha0
        rhoS0 = eos.referenceDensity
        self.rho0 = rhoS0/alpha0
        PS0 = eos.pressure(rhoS0, eps0)
        self.P0 = PS0/alpha0
        if cS0 is None:
            cS0 = eos.soundSpeed(rhosS0, eps0)
        if c0 is None:
            c0 = cS0
        if h0 is None:
            self.h0 = 2.0/nPoints
        else:
            self.h0 = h0
        self.nPoints = nPoints
        self.pistonFrame = pistonFrame
        self.crushCurve = PalphaCrushCurve(rhoS0, self.P0, alpha0, alphat, Pe, Pt, Ps, cS0, c0, n1, n2)
        return

    # Find the shock and elastic wave properties as a function of time
    def waveProperties(self, t):

        # Compute the shock jump conditions
        us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae = computeHugoniotWithPorosity(self.eos, self.rho0, self.eps0, abs(self.vpiston), self.crushCurve)
        ce = self.crushCurve.ce(self.alpha0)
        xs = 1.0 - us*t    # position of the shock
        xe = 1.0 - ce*t    # position of the elastic wave
        # print("Elastic front at ", xe, ", shock front at ", xs)

        # Conditions behind shock
        v1 = -self.vpiston
        h1 = self.h0 * self.rho0*safeInv(rhos)

        # Conditions behind the elastic wave
        v2 = -ue
        h2 = self.h0 * self.rho0*safeInv(rhoe)
        
        return us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae, xs, xe, v1, h1, v2, h2

    # Return the solution profiles as x, v, eps, rho, P, h
    def solution(self, t,
                 x = None):

        # Get the current regional properties
        us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae, xs, xe, v1, h1, v2, h2 = self.waveProperties(t)

        # Did the user specify the x positions?
        xpiston = 1.0 - self.vpiston*t
        if x is None:
            x = np.linspace(-1.0, xpiston, self.nPoints)

        # Generate the solution profiles
        n = len(x)
        v, eps, rho, P, h = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
        for i in range(n):
            if x[i] >= xs:
                v[i], eps[i], rho[i], P[i], h[i] = v1, epss, rhos, Ps, h1
            elif x[i] >= xe:
                v[i], eps[i], rho[i], P[i], h[i] = v2, epse, rhoe, Pe, h2
            else:
                v[i], eps[i], rho[i], P[i], h[i] = 0.0, self.eps0, self.rho0, self.P0, self.h0

        # In the piston frame we actually look like the fluid is flowing into a stationary wall
        if self.pistonFrame:
            x += self.vpiston*t
            v += self.vpiston

        return x, v, eps, rho, P, h

    # Return the alpha solution
    def alpha_solution(self, t,
                       x = None):
        
        # Get the current regional properties
        us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae, xs, xe, v1, h1, v2, h2 = self.waveProperties(t)

        # Did the user specify the x positions?
        if x is None:
            xpiston = 1.0 - self.vpiston*t
            x = np.linspace(-1.0, xpiston, self.nPoints)

        # Generate the solution profiles
        n = len(x)
        alpha = np.zeros(n)
        for i in range(n):
            if x[i] >= xs:
                alpha[i] = alphas
            elif x[i] >= xe:
                alpha[i] = alphae
            else:
                alpha[i] = self.alpha0

        # In the piston frame we actually look like the fluid is flowing into a stationary wall
        if self.pistonFrame:
            x += self.vpiston*t

        return x, alpha

    # Return the sound speed solution
    def soundSpeed_solution(self, t,
                            x = None):
        
        # Get the current regional properties
        us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae, xs, xe, v1, h1, v2, h2 = self.waveProperties(t)
        def _soundSpeed(rhoi, epsi, alphai):
            cS0 = self.eos.soundSpeed(alphai*rhoi, epsi)
            if self.alpha0 > 1.0:
                return cS0 + (alphai - 1.0)/(self.alpha0 - 1.0)*(self.crushCurve.c0 - cS0)
            else:
                return cS0
        c_s = _soundSpeed(rhos, epss, alphas)
        c_e = _soundSpeed(rhoe, epse, alphae)
        c_0 = self.crushCurve.c0

        # Did the user specify the x positions?
        if x is None:
            xpiston = 1.0 - self.vpiston*t
            x = np.linspace(-1.0, xpiston, self.nPoints)

        # Generate the solution profiles
        n = len(x)
        cs = np.zeros(n)
        for i in range(n):
            if x[i] >= xs:
                cs[i] = c_s
            elif x[i] >= xe:
                cs[i] = c_e
            else:
                cs[i] = c_0

        # In the piston frame we actually look like the fluid is flowing into a stationary wall
        if self.pistonFrame:
            x += self.vpiston*t

        return x, cs

#-------------------------------------------------------------------------------
# Main if run directly
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    from Spheral1d import *
    def plotIt(x, y, ylabel, xlabel=r"$x$", ylog=False, title=None, style="k-"):
        fig = newFigure()
        if ylog:
            fig.semilogy(x, y, style)
        else:
            fig.plot(x, y, style)
        fig.set_xlabel(xlabel)
        fig.set_ylabel(ylabel)
        if title:
            fig.set_title(title)
        else:
            fig.set_title(ylabel)
        return

    units = CGuS()
    cgs = CGS()
    mCGSconv = cgs.unitMassKg/units.unitMassKg
    lCGSconv = cgs.unitLengthMeters/units.unitLengthMeters
    tCGSconv = cgs.unitTimeSec/units.unitTimeSec
    PCGSconv = mCGSconv/(tCGSconv*tCGSconv*lCGSconv)
    vCGSconv = lCGSconv/tCGSconv

    eos = TillotsonEquationOfState("aluminum melosh89", etamin = 0.1, etamax = 10.0, units=units)
    #eos = ANEOS("dunite", constants=CGuS())

    alpha0 = 1.275
    rhoS0 = eos.referenceDensity
    eps0 = 0.0
    PS0 = eos.pressure(rhoS0, eps0)
    Pe = 8e8 * PCGSconv
    Ps = 7e9 * PCGSconv
    cS0 = 5.35e5 * vCGSconv
    c0 = 4.11e5 *vCGSconv
    n1 = 0.0
    n2 = 2.0
    alphat = None # (alpha0 - 1.0)*((Ps - Pe)/(Ps - PS0))**2 + 1.0

    alpha_curve = PalphaCrushCurve(rhoS0, PS0/alpha0, alpha0, alphat, Pe, Pe, Ps, cS0, c0, n1, n2)

    P = np.linspace(PS0, 1.5*Ps, 1000)
    alpha = np.array([alpha_curve(x) for x in P])
    from matplotlib import pyplot as plt
    from SpheralMatplotlib import *
    plotIt(P, alpha, xlabel=r"$P$", ylabel=r"$\alpha$", title=r"$\alpha(P)$ crush curve")

    vpiston = -45.8e-3
    #print("Hugoniot solution: ", computeHugoniotWithPorosity(eos, rhoS0/alpha0, eps0, abs(vpiston), alpha_curve))

    solution = PlanarCompactionSolution(eos,
                                        vpiston,
                                        eps0,
                                        alpha0,
                                        alphat,
                                        Pe,
                                        Pe,
                                        Ps,
                                        n1,
                                        n2,
                                        cS0,
                                        c0)
    x, v, eps, rho, P, h = solution.solution(t=3.5)
    xx, alpha = solution.alpha_solution(t=3.5)
    plotIt(x, v, r"$v$")
    plotIt(x, eps, r"$\varepsilon$")
    plotIt(x, rho, r"$\rho$")
    plotIt(x, P, r"$P$", True)
    plotIt(x, alpha, r"$\alpha$")

    # Plot the actual Hugoniot curves
    up = np.geomspace(1e-6, abs(vpiston), 100)
    us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae = computeHugoniotWithPorosity(eos, rhoS0/alpha0, eps0, up, alpha_curve)
    plotIt(rhos, Ps, xlabel=r"$\rho_s$ (g/cc)", ylabel=r"$P_s$ (Mbar)", title=r"$\rho$-$P$ shock Hugoniot", style="ro-")
    plotIt(us, Ps, xlabel=r"$v_s$ (cm/$\mu$sec)", ylabel=r"$P_s$ (Mbar)", title=r"$v_s$-$P$ shock Hugoniot", style="ro-")
    plotIt(up, us, xlabel=r"$v_p$ (cm/$\mu$sec)", ylabel=r"$v_s$ (cm/$\mu$sec)", title=r"$v_p$-$v_s$ shock Hugoniot", style="ro-")
