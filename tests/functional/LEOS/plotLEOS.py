from Spheral3d import *
from SpheralMatplotlib import plotSurface
import matplotlib.pyplot as plt
import numpy as np

# We'll work in CGuS units.
units = PhysicalConstants(0.01,     # Unit length in meters
                          0.001,    # Unit mass in kg
                          1.0e-6)   # Unit time in sec

#-------------------------------------------------------------------------------
# Build an ANEOS SiO2 like thing.
#-------------------------------------------------------------------------------
etaMin, etaMax = 0.2, 5.0
Tmin, Tmax = 1.0e0, 1.0e8
rhoMin, rhoMax = 0.1, 20.0

eosSiO2 =       ANEOS("quartz",          # Material number (offset sequentially from ANEOS.INPUT)
                      constants = units)

# eosTillotsonBasalt = TillotsonEquationOfState(materialName = "basalt",
#                                               etamin = etaMin,
#                                               etamax = etaMax,
#                                               units = units)

eosLEOS_SiO2 = LEOS(materialNumber = 2210,
                    constants = units)

# eosLEOS_Basalt = LEOS(materialNumber = 4358,
#                      constants = units,
#                      dbname = "/usr/workspace/wciuser/PlanetaryDefense/EOS/leos4358.h5")
#                       #leosFileFormat = "HDF")

# eosLEOS_SiO2 = LEOS(materialNumber = 999,
#                    constants = units,
#                    dbname = "eosANEOSSiO2_testLEOS.ascii",
#                    leosFileFormat = "ascii")

# eosLEOS_Basalt = LEOS(materialNumber = 999,
#                       constants = units,
#                       dbname = "./eosBasalt_testLEOS.ascii",
#                       leosFileFormat = "ascii")

#-------------------------------------------------------------------------------
# Plot the pressure, entropy, & sound speed as a function of (rho, eps)
#-------------------------------------------------------------------------------
n = 50
rho = np.geomspace(etaMin * eosSiO2.referenceDensity,
                   etaMax * eosSiO2.referenceDensity,
                   num = n)

stuff = []

for eos, label in (
        (eosSiO2, "ANEOS SiO2"),
        #(eosTillotsonBasalt, "Tillotson"),
        (eosLEOS_SiO2, "LEOS SiO2"),
        #(eosLEOS_Basalt, "LEOS Basalt"),
):
    rho0 = eos.referenceDensity
    epsMin = eos.specificThermalEnergy(rho0, Tmin)
    epsMax = eos.specificThermalEnergy(rho0, Tmax)
    epsOff = min(2.0*epsMin, 0.0)
    eps = np.geomspace(epsMin - epsOff, epsMax, num = n) + epsOff
    T = np.geomspace(Tmin, Tmax, num = n)

    rho_grid, eps_grid = np.meshgrid(rho, eps)
    rho_grid, T_grid = np.meshgrid(rho, T)
    shape = rho_grid.shape
    PA_grid, csA_grid = np.zeros(shape), np.zeros(shape)
    dPdUA_grid, dPdRA_grid = np.zeros(shape), np.zeros(shape)
    sA_grid, gA_grid = np.zeros(shape), np.zeros(shape)
    Teps_grid, epsT_grid, epsTratio_grid = np.zeros(shape), np.zeros(shape), np.zeros(shape)
    K_grid = np.zeros(shape)

    for j in range(n):
        for i in range(n):
            rhoi = rho_grid[j][i]
            epsi = eps_grid[j][i]
            Pi, dPdUi, dPdRi = eos.pressureAndDerivs(rhoi,epsi)
            Pi = eos.pressure(rhoi,epsi)
            PA_grid[j][i] = Pi
            dPdUA_grid[j][i] = dPdUi
            dPdRA_grid[j][i] = dPdRi
            csA_grid[j][i] = eos.soundSpeed(rhoi,epsi)
            sA_grid[j][i] = eos.entropy(rhoi,epsi)
            gA_grid[j][i] = eos.gamma(rhoi,epsi)
            Teps_grid[j][i] = eos.temperature(rhoi,epsi)
            K_grid[j][i] = eos.bulkModulus(rhoi,epsi)

            Ti = T_grid[j][i]
            epsi = eos.specificThermalEnergy(rhoi, Ti)
            Tii = eos.temperature(rhoi, epsi)
            epsT_grid[j][i] = epsi
            epsTratio_grid[j][i] = Tii/Ti

    print("Pressure range for %s    : [%g, %g]" % (label, np.min(PA_grid), np.max(PA_grid)))
    print("dPdU range for %s        : [%g, %g]" % (label, np.min(dPdUA_grid), np.max(dPdUA_grid)))
    print("dPdRho range for %s      : [%g, %g]" % (label, np.min(dPdRA_grid), np.max(dPdRA_grid)))
    print("Sound speed range for %s : [%g, %g]" % (label, np.min(csA_grid), np.max(csA_grid)))
    print("Entropy range for %s     : [%g, %g]" % (label, np.min(sA_grid), np.max(sA_grid)))
    print("Gamma range for %s       : [%g, %g]" % (label, np.min(gA_grid), np.max(gA_grid)))
    print("eps lookup range for %s  : [%g, %g]" % (label, np.min(epsT_grid), np.max(epsT_grid)))
    print("K lookup range for %s    : [%g, %g]" % (label, np.min(K_grid), np.max(K_grid)))
    print("T lookup range for %s    : [%g, %g]" % (label, np.min(Teps_grid), np.max(Teps_grid)))
    print("T(rho,eps)/T range for %s: [%g, %g]" % (label, np.min(epsTratio_grid), np.max(epsTratio_grid)))

    stuff.append(plotSurface(rho_grid, eps_grid - epsOff, PA_grid,
                             xlabel = r"$\rho$ (g/cm$^3$)",
                             ylabel = r"$\varepsilon$ (Mb cm$^2$/g)",
                             zlabel = r"$P$ (dyne)",
                             title = r"Pressure %s" % label))

    stuff.append(plotSurface(rho_grid, eps_grid - epsOff, dPdUA_grid,
                             xlabel = r"$\rho$ (g/cm$^3$)",
                             ylabel = r"$\varepsilon$ (Mb cm$^2$/g)",
                             zlabel = r"$\partial P/\partial \varepsilon$",
                             title = r"$\partial P/\partial \varepsilon$ %s" % label))

    stuff.append(plotSurface(rho_grid, eps_grid - epsOff, dPdRA_grid,
                             xlabel = r"$\rho$ (g/cm$^3$)",
                             ylabel = r"$\varepsilon$ (Mb cm$^2$/g)",
                             zlabel = r"$\partial P/\partial \rho$",
                             title = r"$\partial P/\partial \rho$ %s" % label))

    stuff.append(plotSurface(rho_grid, eps_grid - epsOff, csA_grid,
                             xlabel = r"$\rho$ (g/cm$^3$)",
                             ylabel = r"$\varepsilon$ (Mb cm$^2$/g)",
                             zlabel = r"$c_s$ (cm/sec)",
                             title = r"Sound speed %s" % label))

    stuff.append(plotSurface(rho_grid, eps_grid - epsOff, sA_grid,
                             xlabel = r"$\rho$ (g/cm$^3$)",
                             ylabel = r"$\varepsilon$ (Mb cm$^2$/g)",
                             zlabel = r"$s$",
                             title = r"entropy %s" % label))

    stuff.append(plotSurface(rho_grid, eps_grid - epsOff, gA_grid,
                             xlabel = r"$\rho$ (g/cm$^3$)",
                             ylabel = r"$\varepsilon$ (Mb cm$^2$/g)",
                             zlabel = r"$\gamma$",
                             title = r"gamma %s" % label))

    stuff.append(plotSurface(rho_grid, eps_grid - epsOff, K_grid,
                             xlabel = r"$\rho$ (g/cm$^3$)",
                             ylabel = r"$\varepsilon$ (Mb cm$^2$/g)",
                             zlabel = r"Bulk modulus",
                             title = r"$K$ (dyne) %s" % label))

    stuff.append(plotSurface(rho_grid, eps_grid - epsOff, epsT_grid,
                             xlabel = r"$\rho$ (g/cm$^3$)",
                             ylabel = r"$\varepsilon$ (Mb cm$^2$/g)",
                             zlabel = r"$\varepsilon(\rho, T)$",
                             title = r"eps(rho,T) %s" % label))

    stuff.append(plotSurface(rho_grid, eps_grid - epsOff, Teps_grid,
                             xlabel = r"$\rho$ (g/cm$^3$)",
                             ylabel = r"$\varepsilon$ (Mb cm$^2$/g)",
                             zlabel = r"$T(\rho, \varepsilon)$",
                             title = r"T(rho,eps) %s" % label))

    stuff.append(plotSurface(rho_grid, eps_grid - epsOff, epsTratio_grid,
                             xlabel = r"$\rho$ (g/cm$^3$)",
                             ylabel = r"$\varepsilon$ (Mb cm$^2$/g)",
                             zlabel = r"$T(\rho, \varepsilon/T)$",
                             title = r"T(rho,eps(rho,T))/T %s" % label))
