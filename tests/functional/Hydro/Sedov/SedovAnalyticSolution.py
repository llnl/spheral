#-------------------------------------------------------------------------------
# Provide the analytic solution to the Sedov blastwave problems.
#
# Sedov, L.I. 1959, "Similarity and Dimensional Methods in Mechanics", 210-233.
#-------------------------------------------------------------------------------
from math import *
import numpy as np
import bisect

def sgn(x):
    if x < 0:
        return -1
    else:
        return 1

#-------------------------------------------------------------------------------
# Numerical integration routines, based on Numerical recipes.
#-------------------------------------------------------------------------------
def trapzd(func, a, b, s, n):
    assert b >= a
    assert n >= 1
    if n == 1:
        return 0.5*(b - a)*(func(a) + func(b))
    else:
        it = 2**(n - 2)
        delta = (b - a)/it
        result = 0.0
        for j in range(it):
            result += func(a + (j + 0.5)*delta)
        return 0.5*(s + (b - a)*result/it)

def polint(xa, ya, x):
    n = len(xa)
    assert len(ya) == n
    c = ya[:]
    d = ya[:]
    ns = 0
    dif = abs(x - xa[0])
    for i in range(1, n):
        dift = abs(x - xa[i])
        if dift < dif:
            ns = i
            dif = dift
    y = ya[ns]
    ns -= 1
    for m in range(1, n - 1):
        for i in range(n - m):
            ho = xa[i] - x
            hp = xa[i + m] - x
            w = c[i + 1] - d[i]
            den = ho - hp
            if den == 0.0:
                raise "Failure in polint"
            den = w/den
            d[i] = hp*den
            c[i] = ho*den
        if 2*ns < n - m - 1:
            dy = c[ns + 1]
        else:
            dy = d[ns]
            ns -= 1
        y += dy
    return y, dy

def qromb(func, a, b,
          eps = 1.0e-10,
          maxIters = 50,
          K = 5):
    KM = K - 1
    h = [0.0]*(maxIters + 1)
    s = [0.0]*(maxIters + 1)
    h[0] = 1.0
    for j in range(maxIters):
        jj = j + 1
        s[j] = trapzd(func, a, b, s[j], jj)
        if jj >= K:
            ss, dss = polint(h[j - KM:jj], s[j - KM:jj], 0.0)
            if abs(dss) <= eps*abs(ss):
                return ss
        s[jj] = s[j]
        h[jj] = 0.25*h[j]    
    raise "Too many iterations in qromb"

#-------------------------------------------------------------------------------
# SedovSolution : Main class
#-------------------------------------------------------------------------------
class SedovSolution:

    #---------------------------------------------------------------------------
    # Constructor
    #---------------------------------------------------------------------------
    def __init__(self,
                 nDim,
                 gamma = 1.4,
                 rho0 = 1.0,
                 E0 = 1.0,
                 h0 = 2.01*0.02,
                 nbins = 10001,
                 accuracy = 1e-6):
        self.nu = nDim
        self.gamma = gamma
        self.rho0 = rho0
        self.E0 = E0
        self.h0 = h0
        self.nbins = nbins
        self.accuracy = accuracy

        # Store the alpha exponent constants.
        self.a1 = self.computeAlpha1(self.nu, self.gamma)
        self.a2 = self.computeAlpha2(self.nu, self.gamma)
        self.a3 = self.computeAlpha3(self.nu, self.gamma)
        self.a4 = self.computeAlpha4(self.nu, self.gamma)
        self.a5 = self.computeAlpha5(self.nu, self.gamma)
        self.a6 = self.computeAlpha6(self.nu, self.gamma)
        self.a7 = self.computeAlpha7(self.nu, self.gamma)

        # Compute an internal table with the dimensionless solution.
        self._lam, self._vlam, self._rholam, self._Plam = \
                   self.computeDimensionlessTable(self.nu, self.gamma, nbins)

        # Compute the energy alpha constant.
        self.alpha = self.computeAlpha(nbins)
        assert self.alpha > 0.0

        # Compute E
        self.E = self.E0/self.alpha

        return

    #---------------------------------------------------------------------------
    # Properties at the shock front.
    #---------------------------------------------------------------------------
    def shockState(self, t):

        gamma = self.gamma
        nu = self.nu
        E = self.E
        rho1 = self.rho0

        nu1 = 1.0/(nu + 2.0)
        nu2 = 2.0*nu1
        t1 = t**nu1
        t2 = t**nu2

        vs = nu2 * (E/rho1)**nu1
        if t != 0.0:
            vs /= t**(nu/(nu + 2.0))
        r2 = (E/rho1)**nu1 * t2
        v2 = 2.0/(gamma + 1.0)*vs
        rho2 = (gamma + 1.0)/(gamma - 1.0)*rho1
        P2 = 2.0/(gamma + 1)*rho1*vs*vs

        return vs, r2, v2, rho2, P2

    #---------------------------------------------------------------------------
    # Compute the solution at the given positions.
    #---------------------------------------------------------------------------
    def solution(self, t,
                 r = None):

        gamma = self.gamma
        gam1 = gamma - 1.0
        nu = self.nu

        vs, r2, v2, rho2, P2 = self.shockState(t)

        if r is None:
            r = np.linspace(0.0, r2, 101, endpoint=True)
            r = np.append(r, 2.0*r2)
            #r = np.array([0.01*r2*i for i in range(101)])

        v, rho, P, u, h, A = np.zeros_like(r), np.zeros_like(r), np.zeros_like(r), np.zeros_like(r), np.zeros_like(r), np.zeros_like(r)

        for i, ri in enumerate(r):
            if abs(ri) < r2:
                vi, rhoi, Pi = self.lookupSolution(abs(ri/r2))
            else:
                vi, rhoi, Pi = 0.0, self.rho0/rho2, 0.0
            v[i] = vi*v2*sgn(ri)
            rho[i] = rhoi*rho2
            P[i] = Pi*P2
            u[i] = P[i]/(gam1*rho[i] + 1.0e-50)
            h[i] = self.h0 * self.rho0/(rho[i] + 1.0e-50)
            A[i] = P[i]/max(1.0e-30, rho[i])**gamma

        return r, v, u, rho, P, A, h

    #---------------------------------------------------------------------------
    # Compute the derivatives
    # Numerically estimates and returns (dvdt, dudt, drhodt, divv)
    #---------------------------------------------------------------------------
    def derivatives(self, t,
                    r = None,
                    rmin = None,
                    rmax = None,
                    numSolution = 1000):
        
        gamma = self.gamma
        gam1 = gamma - 1.0
        nu = self.nu

        vs, r2, v2, rho2, P2 = self.shockState(t)

        if r is None:
            r = np.array([0.01*r2*i for i in range(101)])

        if rmin is None:
            rmin = min(r)
        if rmax is None:
            rmax = max(r)

        # Make the solution grid in r to compute derivatives on
        rsol = np.linspace(rmin, rmax, numSolution, endpoint=True)

        # Compute the standard solution on the fine sampled grid
        rdummy, v, u, rho, P, A, h = self.solution(t, rsol)

        gradP = np.gradient(P, rsol)
        divv = np.gradient(v, rsol)
        thpt = np.zeros_like(divv)
        if self.nu == 2:
            divv += np.divide(v, rsol, where = (rsol != 0.0))
        elif self.nu == 3:
            divv += np.divide(v, rsol*rsol, where = (rsol != 0.0))

        dvdt = -np.divide(gradP, rho, where = (rho != 0.0))
        dudt = -np.divide(P, rho, where = (rho != 0.0)) * divv
        drhodt = -rho*divv

        return (np.interp(r, rsol, dvdt),
                np.interp(r, rsol, dudt),
                np.interp(r, rsol, drhodt),
                np.interp(r, rsol, divv))

    #---------------------------------------------------------------------------
    # alpha1
    #---------------------------------------------------------------------------
    def computeAlpha1(self, nu, gamma):
        return (((nu + 2.0)*gamma)/(2.0 + nu*(gamma - 1.0)) *
                (2.0*nu*(2.0 - gamma)/(gamma*(nu + 2)**2) - self.computeAlpha2(nu, gamma)))

    #---------------------------------------------------------------------------
    # alpha2
    #---------------------------------------------------------------------------
    def computeAlpha2(self, nu, gamma):
        return (1.0 - gamma)/(2.0*(gamma - 1.0) + nu)

    #---------------------------------------------------------------------------
    # alpha3
    #---------------------------------------------------------------------------
    def computeAlpha3(self, nu, gamma):
        return nu/(2.0*(gamma - 1.0) + nu)

    #---------------------------------------------------------------------------
    # alpha4
    #---------------------------------------------------------------------------
    def computeAlpha4(self, nu, gamma):
        return self.computeAlpha1(nu, gamma)*(nu + 2.0)/(2.0 - gamma)

    #---------------------------------------------------------------------------
    # alpha5
    #---------------------------------------------------------------------------
    def computeAlpha5(self, nu, gamma):
        return 2.0/(gamma - 2.0)

    #---------------------------------------------------------------------------
    # alpha6
    #---------------------------------------------------------------------------
    def computeAlpha6(self, nu, gamma):
        return gamma/(2.0*(gamma - 1.0) + nu)

    #---------------------------------------------------------------------------
    # alpha7
    #---------------------------------------------------------------------------
    def computeAlpha7(self, nu, gamma):
        return (2.0 + nu*(gamma - 1))*self.computeAlpha1(nu, gamma)/(nu*(2.0 - gamma))

    #---------------------------------------------------------------------------
    # lam (r/r2) 
    #---------------------------------------------------------------------------
    def lam(self, V):
        nu = self.nu
        gamma = self.gamma
        a1 = self.a1
        a2 = self.a2

        return ((0.25*(nu + 2.0)*(gamma + 1)*V)**(-2.0/(2.0 + nu)) *
                ((gamma + 1)/(gamma - 1.0) * (0.5*(nu + 2.0)*gamma*V - 1.0))**(-a2) *
                ((nu + 2.0)*(gamma + 1.0)/((nu + 2.0)*(gamma + 1.0) - 2.0*(2.0 + nu*(gamma - 1.0))) *
                 (1.0 - 0.5*(2.0 + nu*(gamma - 1.0))*V))**(-a1)
                )

    #---------------------------------------------------------------------------
    # vlambda (v/v2)
    #---------------------------------------------------------------------------
    def vlambda(self, V):
        nu = self.nu
        gamma = self.gamma

        return 0.25*(nu + 2.0)*(gamma + 1.0) * V * self.lam(V)

    #---------------------------------------------------------------------------
    # rholambda (rho/rho2)
    #---------------------------------------------------------------------------
    def rholambda(self, V):
        nu = self.nu
        gamma = self.gamma
        a3 = self.a3
        a4 = self.a4
        a5 = self.a5

        return (((gamma + 1.0)/(gamma - 1.0)*(0.5*(nu + 2.0)*gamma*V - 1.0))**a3 *
                ((gamma + 1.0)/(gamma - 1.0)*(1.0 - 0.5*(nu + 2.0)*V))**a5 *
                ((nu + 2.0)*(gamma + 1.0)/((2.0 + nu)*(gamma + 1.0) - 2.0*(2.0 + nu*(gamma - 1.0))) *
                 (1.0 - 0.5*(2.0 + nu*(gamma - 1.0))*V))**a4)

    #---------------------------------------------------------------------------
    # Plambda (P/P2)
    #---------------------------------------------------------------------------
    def Plambda(self, V):
        nu = self.nu
        gamma = self.gamma
        a1 = self.a1
        a4 = self.a4
        a5 = self.a5

        return ((0.25*(nu + 2.0)*(gamma + 1.0)*V)**(2.0*nu/(2.0 + nu)) *
                ((gamma + 1.0)/(gamma - 1.0)*(1.0 - 0.5*(nu + 2.0)*V))**(a5 + 1.0) *
                ((nu + 2.0)*(gamma + 1.0)/((nu + 2.0)*(gamma + 1.0) - 2.0*(2.0 + nu*(gamma - 1.0))) *
                 (1.0 - 0.5*(2.0 + nu*(gamma - 1.0))*V))**(a4 - 2.0*a1))

    #---------------------------------------------------------------------------
    # The range of the dimensionless velocity variable V.
    #---------------------------------------------------------------------------
    def Vrange(self, nu, gamma):
        assert gamma > 1.0
        assert nu in (1, 2, 3)
        if (nu in (1, 2)) or (nu == 3 and gamma < 7.0):
            return (2.0/((nu + 2.0)*gamma), 4.0/((nu + 2.0)*(gamma + 1.0)))
        else:
            return (4.0/(5.0*(gamma + 1.0)), 2.0/5.0)

    #---------------------------------------------------------------------------
    # The dimension dependent volume constant.
    #---------------------------------------------------------------------------
    def Anu(self, nu):
        if nu == 1:
            return 2.0
        elif nu == 2:
            return 2.0*pi
        elif nu == 3:
            return 4.0*pi
        else:
            assert False

    #---------------------------------------------------------------------------
    # Compute a tabular form of the dimensionless solution.
    #---------------------------------------------------------------------------
    def computeDimensionlessTable(self, nu, gamma, nbins):
        assert nbins > 1

        gam1 = gamma - 1.0

        Vmin, Vmax = self.Vrange(nu, gamma)
        dV = (Vmax - Vmin)/(nbins - 1)

        lam = []
        vlam = []
        rholam = []
        Plam = []
        for i in range(nbins):
            V = Vmin + i*dV
            lam.append(self.lam(V))
            vlam.append(self.vlambda(V))
            rholam.append(self.rholambda(V))
            Plam.append(self.Plambda(V))

        assert len(lam) == nbins
        assert len(vlam) == nbins
        assert len(rholam) == nbins
        assert len(Plam) == nbins

        assert min([x >= 0.0 and x <= 1.0 for x in lam]) == True
        assert lam[0] < 1.0e-10 and lam[-1] > 1.0 - 1.0e-10
        lam[0] = 0.0
        lam[-1] = 1.0

        return lam, vlam, rholam, Plam

    #---------------------------------------------------------------------------
    # Peform a linear interpolation into the given table.
    #---------------------------------------------------------------------------
    def interp(self, x, xtab, ytab, imin, imax):
        return ytab[imin] + ((ytab[imax] - ytab[imin])/
                             (xtab[imax] - xtab[imin] + 1.0e-50)*
                             (x - xtab[imin]))

    #---------------------------------------------------------------------------
    # Interpolate for the solution at a given lambda \in [0,1].
    #---------------------------------------------------------------------------
    def lookupSolution(self, x):
        assert x >= 0.0 and x <= 1.0

        # Bracket this x in the table.
        imin = max(0, min(self.nbins - 1, bisect.bisect(self._lam, x) - 1))
        imax = min(imin + 1, self.nbins - 1)
        assert imin >= 0 and imin < self.nbins
        assert imax >= 0 and imax < self.nbins

        # Now we can interpolate the values.
        v = self.interp(x, self._lam, self._vlam, imin, imax)
        rho = self.interp(x, self._lam, self._rholam, imin, imax)
        P = self.interp(x, self._lam, self._Plam, imin, imax)

        return v, rho, P

    #---------------------------------------------------------------------------
    # Numerically integrate the alpha constant.
    #---------------------------------------------------------------------------
    def computeAlpha(self, accuracy):
        assert accuracy > 0.0
        gamma = self.gamma
        nu = self.nu
        Anu = self.Anu(nu)
        thpt = qromb(self.func, 0.0, 1.0, accuracy)
        return 8.0*Anu/((gamma - 1.0)*(gamma + 1.0)*(nu + 2.0)**2) * thpt

    #---------------------------------------------------------------------------
    # The integrand function for integrating alpha
    #---------------------------------------------------------------------------
    def func(self, x):
        assert x >= 0.0 and x <= 1.0
        v0, rho0, P0 = self.lookupSolution(x)
        return (rho0*v0**2 + P0) * x**(self.nu - 1)

##    #---------------------------------------------------------------------------
##    # Numerically integrate the alpha constant.
##    #---------------------------------------------------------------------------
##    def computeAlpha(self, nbins, f):
##        assert nbins > 1
##        nu = self.nu
##        nu1 = nu - 1
##        gamma = self.gamma
##        Anu = self.Anu(nu)

##        # First figure out our geometric ratio sizing base.
##        fsum = 0.0
##        for i in xrange(nbins):
##            fsum += f**i
##        assert fsum > 0.0
##        dlambda0 = 1.0/fsum

##        result = 0.0
##        lambdasum = 0.0
##        for i in xrange(nbins):
##            dlambdai = dlambda0 * f**i
##            lam0 = lambdasum
##            lam1 = min(1.0, lambdasum + dlambdai)
##            assert lam0 >= 0.0 and lam0 <= 1.0
##            assert lam1 >= 0.0 and lam1 <= 1.0
##            v0, rho0, P0 = self.lookupSolution(lam0)
##            v1, rho1, P1 = self.lookupSolution(lam0)
##            val0 = (rho0*v0**2 + P0)*lam0**nu1
##            val1 = (rho1*v1**2 + P1)*lam1**nu1
##            result += 0.5*(val0 + val1)*(lam1 - lam0)
##            lambdasum = lam1
##        assert abs(1.0 - lambdasum) < 1.0e-10
##        result *= 8.0*Anu/((gamma - 1.0)*(gamma + 1.0)*(nu + 2.0)**2)

##        return result # 1.1*result

##    #---------------------------------------------------------------------------
##    # Numerically integrate the alpha constant.
##    #---------------------------------------------------------------------------
##    def computeAlpha(self, nbins, f):
##        assert nbins > 1
##        nu = self.nu
##        nu1 = nu - 1
##        gamma = self.gamma
##        Anu = self.Anu(nu)

##        # First figure out our geometric ratio sizing base.
##        fsum = 0.0
##        for i in xrange(nbins):
##            fsum += f**i
##        assert fsum > 0.0
##        dlambda0 = 1.0/fsum

##        result = 0.0
##        lambdasum = 0.0
##        for i in xrange(nbins):
##            dlambdai = dlambda0 * f**i
##            lam0 = lambdasum
##            lam1 = min(1.0, lambdasum + dlambdai)
##            assert lam0 >= 0.0 and lam0 <= 1.0
##            assert lam1 >= 0.0 and lam1 <= 1.0
##            v0, rho0, P0 = self.lookupSolution(lam0)
##            v1, rho1, P1 = self.lookupSolution(lam0)
##            val0 = (rho0*v0**2 + P0)*lam0**nu1
##            val1 = (rho1*v1**2 + P1)*lam1**nu1
##            result += 0.5*(val0 + val1)*(lam1 - lam0)
##            lambdasum = lam1
##        assert abs(1.0 - lambdasum) < 1.0e-10
##        result *= 8.0*Anu/((gamma - 1.0)*(gamma + 1.0)*(nu + 2.0)**2)

##        return result # 1.1*result

##    #---------------------------------------------------------------------------
##    # Numerically integrate the alpha constant.
##    #---------------------------------------------------------------------------
##    def computeAlpha(self, nbins):
##        assert nbins > 1
##        nu = self.nu
##        nu1 = nu - 1
##        gamma = self.gamma
##        Anu = self.Anu(nu)

##        Vmin, Vmax = self.Vrange(nu, gamma)
##        dV = (Vmax - Vmin)/nbins

##        result = 0.0
##        for i in xrange(nbins):
##            V0 = Vmin + i*dV
##            V1 = min(Vmax, Vmin + (i + 1)*dV)
##            assert V0 >= Vmin and V0 <= Vmax
##            assert V1 >= Vmin and V1 <= Vmax
##            lam0 = self.lam(V0)
##            lam1 = self.lam(V1)
##            dlambda = lam1 - lam0
##            assert dlambda > 0.0
##            val0 = (self.rholambda(V0)*self.vlambda(V0)**2 + self.Plambda(V0))*lam0**nu1
##            val1 = (self.rholambda(V1)*self.vlambda(V1)**2 + self.Plambda(V1))*lam1**nu1
##            result += 0.5*(val0 + val1)*dlambda
##        result *= 8.0*Anu/((gamma - 1.0)*(gamma + 1.0)*(nu + 2.0)**2)

##        return result
