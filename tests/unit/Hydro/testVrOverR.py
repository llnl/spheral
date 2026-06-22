#-------------------------------------------------------------------------------
# Test out using analytic integral of vr/r vs single point derivative
#-------------------------------------------------------------------------------
from Spheral import *
import numpy as np
import matplotlib.pyplot as plt
from SpheralMatplotlib import *
from SpheralTestUtilities import *

#-------------------------------------------------------------------------------
# Analytic integral of vr/r over a timestep
#-------------------------------------------------------------------------------
def integrate_vr_over_r(vr,
                        r,
                        ar,
                        dt):
    tiny = 1.0e-6
    if r < 0.0:
        r = -r;
        vr = -vr;
        ar = -ar;

    if (fuzzyEqual(vr, 0.0, tiny) and fuzzyEqual(ar, 0.0, tiny)):
        return 0.0;

    def D(t):
        return r + vr*t + ar*t*t

    def F0(t):
        return 0.5*log(abs(D(t)*safeInvVar(r)))

    def F1(t):
        X = 4.0*ar*r - vr*vr
        X12 = sqrt(abs(X))
        if (X > 0.0):
            result = 2.0/X12*(atan2(vr + 2.0*ar*t, X12) - atan2(vr, X12))
            return result
        elif (X < 0.0):
            if (fuzzyEqual(vr + 2.0*ar*t - X12, 0.0, tiny) or fuzzyEqual(vr + X12, tiny)):
                return 0.0;
            thpt = abs((vr + 2.0*ar*t - X12)*(vr + X12)*safeInvVar((vr + 2.0*ar*t + X12)*(vr - X12), tiny))
            # if (fuzzyEqual(thpt, 0.0, tiny)) return 0.0;
            result = log(thpt)/X12;
            return result
        else:
            assert X == 0.0
            result = 2.0*ar*(safeInvVar(vr) - safeInvVar(vr + 2.0*ar*t))*safeInvVar(ar);
            return result

    result = F0(dt) + 0.5*vr*F1(dt)
    return result/dt

#-------------------------------------------------------------------------------
# Generate a grid of comparison values
#-------------------------------------------------------------------------------
N = 100
rmin, rmax = 1.0/N, 1.0 + 1.0/N
ar0 = -0.1
vr0 = -0.2
dt = 1e-3
r = np.linspace(rmin, rmax, N, endpoint=True)
ar = np.full(r.shape, ar0)
vr = np.full(r.shape, vr0)
vr_over_r0 = vr/r
vr_over_r1 = np.array([SPHRZ.integrate_vr_over_r(vr[i], r[i], ar[i], dt) for i in range(N)])

fig0 = newFigure()
fig0.plot(r, vr_over_r0, "k-", label="$v_r/r$")
fig0.plot(r, vr_over_r1, "b-.", label="$\\int v_r/r$")
fig0.set_xlabel("$r$")
fig0.set_ylabel("$v_r/r$")
fig0.axes.legend()

fig1 = newFigure()
fig1.plot(r, (vr_over_r1 - vr_over_r0)/vr_over_r0, "k-", label="$(\\int v_r/r dr - v_r/r)/(v_r/r)$")
fig1.set_xlabel("$r$")
fig1.axes.legend()

fig2 = newFigure()
fig2.plot(r, vr_over_r1/vr_over_r0, "k-", label="$(\\int v_r/r dr)/(v_r/r)$")
fig2.set_xlabel("$r$")
fig2.axes.legend()
