#-------------------------------------------------------------------------------
# Generalized version of computeHugoniot to include porosity, which can lead to
# a new elastic porosity compression state in front of the shock front.
#
# Based on the discussion in the section on Hugoniots and Isentropes from
# Herrmann, W. (1969). Constitutive Equation for the Dynamic Compaction of 
# Ductile Porous Materials. Journal of Applied Physics, 40(6), 2490–2499.
#
# Also see
# Jutzi, M., Benz, W., & Michel, P. (2008). Numerical simulations of impacts
# involving porous bodies.I. Implementing sub-resolution porosity in a 3D SPH
# hydrocode. Icarus, 198(1), 242–255.
#
# Compute Hugoniot shock jump values for the given equation of state based on an
# initial (density, energy, alpha0, upiston) state.  upiston can be a single value or
# a numpy array.
# Arguments:
#       eos: the equation of state to profile
#      rho0: initial mass density (bulk or porous value)
#      eps0: initial specific thermal energy
#   upiston: piston velocity (single value or numpy array)
# 
# Returns:
#        us: shock velocity
#      rhos: post-shock density
#      epss: post-shock specific thermal energy
#        Ps: post-shock pressure
# 
# If upiston is a numpy array then the returned values are also arrays corresponding
# to each piston value.
#-------------------------------------------------------------------------------
import mpi
import numpy as np
from numpy import linalg as LA
import scipy.integrate, scipy.optimize, scipy.stats
from computeHugoniot import solve

#-------------------------------------------------------------------------------
# Couple of useful functions
#-------------------------------------------------------------------------------
def sgn(x):
    return -1.0 if x < 0.0 else 1.0

def safeInv(x, fuzz = 1e-50):
    return sgn(x)/max(fuzz, abs(x))

#-------------------------------------------------------------------------------
# The P-alpha crush curve from Jutzi 2008
#-------------------------------------------------------------------------------
class PalphaCrushCurve:

    def __init__(self,
                 rhoS0,
                 P0,
                 alpha0,   # alpha(P0)
                 alphat,
                 Pe,
                 Pt,
                 Ps,
                 cS0,
                 c0,
                 n1,
                 n2):
        self.rhoS0 = rhoS0
        self.P0 = P0
        self.alpha0 = alpha0
        self.alphat = alphat
        self.Pe = Pe
        self.Pt = Pt
        self.Ps = Ps
        self.cS0 = cS0
        self.c0 = c0
        self.K0 = cS0*cS0*rhoS0
        self.n1 = n1
        self.n2 = n2

        # Find alphae = alpha(Pe)
        self.alphae = alpha0   # Starting point
        if alpha0 > 1.0:
            last_alphae = 0.0
            iter = 0
            while abs(self.alphae - last_alphae) > 1.0e-15 and iter < 1000:
                iter += 1
                last_alphae = self.alphae
                self.alphae = scipy.integrate.solve_ivp(self.Dalpha_elasticDP,
                                                        t_span = [self.P0, self.Pe],
                                                        y0 = [self.alpha0],
                                                        t_eval = [self.Pe]).y[0][0]

        if self.alphat is None:
            self.alphat = self.alphae    # Reduces to Eq 8 in Jutzi 2008

        assert self.rhoS0 > 0.0
        assert 1.0 <= self.alphat <= self.alphae 
        assert self.Pe <= self.Pt <= self.Ps
        assert self.c0 <= self.cS0
        return

    def h(self, alpha):
        assert self.alphae >= 1.0 and self.c0 <= self.cS0, "alphae={}, c0={}, cS0={}".format(self.alphae, self.c0, self.cS0)
        return 1.0 + (alpha - 1.0)*(self.c0 - self.cS0)/max(1.0e-20, self.cS0*(self.alphae - 1.0))

    def Dalpha_elasticDP(self, P, alpha):
        return alpha*alpha/self.K0*(1.0 - 1.0/self.h(alpha)**2)

    def ce(self, alpha):
        return self.h(alpha)*self.cS0

    def alphaElastic(self, P):
        "Just the elastic curve response"
        stuff = scipy.integrate.solve_ivp(self.Dalpha_elasticDP,
                                          t_span = [self.P0, P],
                                          y0 = [self.alpha0],
                                          t_eval = [P])
        return stuff.y[0][0]

    def __call__(self, P):
        if P <= self.P0:
            return self.alpha0

        if (P >= self.Ps):

            # Solid limit, all porosity squeezed out
            return 1.0
        
        elif P < self.Pe:
            assert self.c0 < self.cS0

            # Elastic regime -- integrate elastic Dalpha/DP equation
            stuff = scipy.integrate.solve_ivp(self.Dalpha_elasticDP,
                                              t_span = [self.P0, P],
                                              y0 = [self.alpha0],
                                              t_eval = [P])
            return stuff.y[0][0]

        else:

            # Plastic limit
            if P < self.Pt:
                return (self.alphae - self.alphat) * ((self.Pt - P)/(self.Pt - self.Pe))**self.n1 + (self.alphat - 1.0) * ((self.Ps - P)/(self.Ps - self.Pe))**self.n2 + 1.0
            else:
                return (self.alphat - 1.0) * ((self.Ps - P)/(self.Ps - self.Pe))**self.n2 + 1.0
                
#-------------------------------------------------------------------------------
# Hugoniot solution with porosity
#-------------------------------------------------------------------------------
def computeHugoniotWithPorosity(eos, rho0, eps0, upiston, crushCurve):

    if type(upiston) is float:
        upiston = np.array([upiston])
    assert np.min(upiston) >= 0.0
    n = len(upiston)

    # How should we query the single pressure response from the EOS?  Some
    # Spheral equations of state support this, but others require the Field
    # interface
    if hasattr(eos, "pressure"):
        Pfunc = eos.pressure
    else:
        raise RuntimeError("computeHugoniotWithPorosity does not work yet with EOS's that can't compute a single pressure")

    # Extract a few constants from the crush curve
    alpha0 = crushCurve.alpha0
    Pe = crushCurve.Pe
    alphae = crushCurve.alphae
    ce = crushCurve.ce(alpha0)    # Elastic wave speed
    P0 = Pfunc(alpha0*rho0, eps0)/alpha0

    # print("Initial conditions: u0: ", 0.0, "\n",
    #       "                  rho0: ", rho0, "\n",
    #       "                  eps0: ", eps0, "\n",
    #       "                    P0: ", P0, "\n",
    #       "                alpha0: ", alpha0, "\n")

    # Functor to help us solve the Rankine-Hugoniot jump relations including a porosity
    class RankineHugoniotJumpRelations:
        def __init__(self, upiston, u0, rho0, eps0, P0, alpha0, alphaPfunc):
            self.upiston = upiston
            self.u0 = u0
            self.rho0 = rho0
            self.eps0 = eps0
            self.P0 = P0
            self.alpha0 = alpha0
            self.alphaPfunc = alphaPfunc
            return

        def __call__(self, args):
            us, eps1, alpha1 = args
            m1 = self.rho0*(us - self.u0)            # mass/time
            rho1 = m1*safeInv(us - self.upiston)
            P1 = Pfunc(alpha1*rho1, eps1)/alpha1
            alpha1_new = self.alphaPfunc(P1)
            return np.array([m1*(self.upiston - self.u0) - (P1 - self.P0),                                               # Conservation of momentum
                             m1*(eps1 - self.eps0 + 0.5*(self.upiston - self.u0)**2) - P1*(self.upiston - self.u0),      # Conservation of energy
                             alpha1_new - alpha1])                                                                       # Convergence of alpha(P)

        def norm(self, args):
            return LA.norm(self(args))
        
        # def jacobian(self, args):
        #     if PDfunc is None:
        #         return None
        #     us, eps1, alpha1 = args
        #     m1 = self.rho0*(us - self.u0)            # mass/time
        #     rho1 = m1*safeInv(us - self.upiston)
        #     PS1, dPSdE1, dPSdR1 = PDfunc(alpha1*rho1, eps1)
        #     P1, dPdE1, dPdR1 = Pfunc(alpha1*rho1, eps1)/alpha1
        #     dUSdR = -self.upiston*self.rho0/max(rho1 - self.rho0)**2
        #     return np.array([self.rho0*(self.upiston - self.u0) - dPdR1*safeInv(dUSdR),    # d(momentum)/dUS
        #                      -dPdE1,                                                       # d(momentum)/dEPS

    # Solve for elastic wave properties (wave moves at speed ce).
    # We need to know what is the maximum piston speed that just generates an elastic wave,
    # but does not create a shock. This is denoted upe.
    #   Known: Pe, ce, alphae
    # Unknown: ue, rhoe, epse, upe
    # Solve for upe
    def elasticLimitConditions(upi):
        m1 = rho0*ce
        rho1 = m1*safeInv(ce - upi)
        eps1 = eps0 - 0.5*upi*upi + (Pe - P0)*upi*safeInv(m1)
        P1 = Pfunc(alphae*rho1, eps1)/alphae
        return P1 - Pe

    # Find the critical piston velocity that just gives us the elastic wave conditions
    stuff = scipy.optimize.root_scalar(elasticLimitConditions, bracket = (0.0, ce))
    upe = stuff.root
    # print("Found elastic limit piston velocity upe: ", upe, "\n", stuff)

    # With upe, recover the full elastic limit conditions
    rhoe = rho0*ce*safeInv(ce - upe)
    epse = eps0 - 0.5*upe*upe + (Pe - P0)*upe/(rho0*ce)
    # print("Elastic conditions:  ce = ", ce, "\n",
    #       "                   rhoe = ", rhoe, "\n",
    #       "                   epse = ", epse, "\n",
    #       "                     Pe = ", Pe, Pfunc(alphae*rhoe, epse)/alphae, abs(Pfunc(alphae*rhoe, epse)/alphae - Pe)/Pe, "\n",
    #       "                 alphae = ", alphae, "\n",
    #       "                    upe = ", upe)

    # Prepare to return the arrays of values.  We return these in the same frame as the piston velocity was given, so presumably lab
    # e => elastic region
    # s => post shock region
    US, RHOS, EPSS, PS, ALPHAS = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    UE, RHOE, EPSE, PE, ALPHAE = np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n)
    for i, up in enumerate(upiston):

        if up <= upe:

            # There is no shock, just the elastic wave
            # print("NO SHOCK -- ELASTIC WAVE ONLY: ", up, upe)
            u1, eps1, alpha1 = solve(func = RankineHugoniotJumpRelations(up, 0.0, rho0, eps0, P0, alpha0, crushCurve.alphaElastic),
                                     initial_guess = (1.5*up, 0.5*up**2, alpha0),
                                     bounds = [(up, ce),         # u
                                               (eps0, epse),     # eps
                                               (1.0, alphae)])   # alpha
            rho1 = rho0*u1*safeInv(u1 - up)
            UE[i] = up
            RHOE[i] = rho1
            EPSE[i] = eps1
            PE[i] = Pfunc(alpha1*rho1, eps1)/alpha1
            ALPHAE[i] = alpha1

        else:

            # Now we have the full elastic wave properties, so solve for the shock jump conditions encountering the
            # post-elastic wave initial conditions.
            us, epss, alphas = solve(func = RankineHugoniotJumpRelations(up, upe, rhoe, epse, Pe, alphae, crushCurve),
                                     initial_guess = (1.5*up,                     # us
                                                      0.5*up**2,                  # epss
                                                      1.0 + 0.5*(alphae - 1.0)),  # alphas
                                     bounds = [(up, 2.0*ce),       # us
                                               (0.0, np.inf),      # epss
                                               (1.0, alpha0)],     # alphas
                                     verbose = False)
            rhos = rhoe*(us - upe)*safeInv(us - up)
            # print("  Shock conditions:  us = ", us, "\n",
            #       "                   rhos = ", rhos, "\n",
            #       "                   epss = ", epss, "\n",
            #       "                 alphas = ", alphas, "\n",
            #       "                     up = ", up)

            # If the shock speed exceeds ce, then the shock front overtakes the elastic wave and there is no
            # elastic region
            if us >= ce:
                # print("NO ELASTIC WAVE -- SHOCK SOLUTION ONLY: ", up, upe, us, ce)
                us, epss, alphas = solve(func = RankineHugoniotJumpRelations(up, 0.0, rho0, eps0, P0, alpha0, crushCurve),
                                         initial_guess = (1.5*up, 0.5*up**2, 1.0),
                                         bounds = [(up, 2.0*ce),       # us
                                                   (0.0, np.inf),      # epss
                                                   (1.0, alpha0)])     # alphas
                rhos = rho0*us*safeInv(us - up)
                # print("  Shock conditions:  us = ", us, "\n",
                #       "                   rhos = ", rhos, "\n",
                #       "                   epss = ", epss, "\n",
                #       "                 alphas = ", alphas)

                US[i] = us
                RHOS[i] = rhos
                EPSS[i] = epss
                PS[i] = Pfunc(alphas*rhos, epss)/alphas
                ALPHAS[i] = alphas

            else:

                # There is both an elastic wave and shocked region
                # print("SHOCK AND ELASTIC REGION: ", up, upe, us, ce)
                US[i] = us
                RHOS[i] = rhos
                EPSS[i] = epss
                PS[i] = Pfunc(alphas*rhos, epss)/alphas
                ALPHAS[i] = alphas
                UE[i] = upe
                RHOE[i] = rhoe
                EPSE[i] = epse
                PE[i] = Pe
                ALPHAE[i] = alphae

    result = US, RHOS, EPSS, PS, ALPHAS, UE, RHOE, EPSE, PE, ALPHAE
    if n == 1:
        return [x[0] for x in result]
    else:
        return result
 
#-------------------------------------------------------------------------------
# Main if run directly
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    from Spheral1d import *
    from matplotlib import pyplot as plt
    from SpheralMatplotlib import *
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

    # Plot the crush curve
    P = np.linspace(PS0, 1.5*Ps, 1000)
    alpha = np.array([alpha_curve(x) for x in P])
    plotIt(P, alpha, xlabel=r"$P$", ylabel=r"$\alpha$", title=r"$\alpha(P)$ crush curve")

    # Plot the Hugoniot curves
    vpiston = -45.8e-3
    up = np.linspace(1e-6, 0.05, 50)
    us, rhos, epss, Ps, alphas, ue, rhoe, epse, Pe, alphae = computeHugoniotWithPorosity(eos, rhoS0/alpha0, eps0, up, alpha_curve)
    plotIt(rhos, Ps, xlabel=r"$\rho_s$ (g/cc)", ylabel=r"$P_s$ (Mbar)", title=r"$\rho$-$P$ shock Hugoniot", style="ro-")
    plotIt(us, Ps, xlabel=r"$v_s$ (cm/$\mu$sec)", ylabel=r"$P_s$ (Mbar)", title=r"$v_s$-$P$ shock Hugoniot", style="ro-")
    plotIt(up, us, xlabel=r"$v_p$ (cm/$\mu$sec)", ylabel=r"$v_s$ (cm/$\mu$sec)", title=r"$v_p$-$v_s$ shock Hugoniot", style="ro-")
