#-------------------------------------------------------------------------------
# Analysis methods to plot the Sedov RZ results
#-------------------------------------------------------------------------------
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import numpy as np
from SedovAnalyticSolution import SedovSolution

plt.rc("font", size=18)

# Parameters
gamma = 5.0/3.0
E0 = 1.0
rho0 = 1.0
nsol = 1000
rshock = 0.8
comp_style, comp_color = "b-", "blue"
noncomp_style, noncomp_color = "r-", "red"

#-------------------------------------------------------------------------------
# Method to read our tabulated data
#-------------------------------------------------------------------------------
def readGnuFile(filename):
    with open(filename, 'r') as f:
        line = f.readline()
        labels = line.split(" ")[2:-1]
        labels = [x.replace("'", "") for x in labels]
        nlabels = len(labels)
        vals = []
        for line in f:
            if line[0] != "#":
                lvals = [eval(_x) for _x in line.split()]
                assert len(lvals) == nlabels
                vals.append(lvals)
    return labels, np.array(vals)

#-------------------------------------------------------------------------------
# Plot averaged values with sigma
#-------------------------------------------------------------------------------
def plotAverage(x, y,
                nbins = 100,
                plot = None,
                lineStyle = "b-",
                markerSize = 4,
                fillcolor = "blue",
                fillalpha = 0.3,
                linewidth = 2,
                lineTitle = None,
                kwords = {}):

    if plot is None:
        plot = plt.figure().add_subplot(111)

    # Use scipy to computed the binned mean and standard deviation
    bin_means, bin_edges, _ = binned_statistic(x, y, statistic='mean', bins=nbins)
    bin_stds, _, _ = binned_statistic(x, y, statistic='std', bins=nbins)

    # Calculate bin centers for plotting
    bin_centers = 0.5*(bin_edges[:-1] + bin_edges[1:])

    # Plot shaded 1-sigma region
    if fillcolor:
        plot.fill_between(bin_centers, bin_means - bin_stds, bin_means + bin_stds, 
                          color=fillcolor, alpha=fillalpha, **kwords)

    # Plot mean line
    if lineStyle:
        plot.plot(bin_centers, bin_means, lineStyle, lw=linewidth, label=lineTitle, **kwords)

    return plot

#-------------------------------------------------------------------------------
# Plot comparison curves
#-------------------------------------------------------------------------------
def plotIt(labels, vals1, vals2,
           xkey, ykey,
           xans, yans):

    result = plt.figure().add_subplot(111)
    ix, iy = labels.index(xkey), labels.index(ykey)
    plotAverage(vals1[:,ix], vals1[:,iy],
                plot = result,
                lineStyle = comp_style,
                fillcolor = comp_color,
                lineTitle = "Compatible")
    plotAverage(vals2[:,ix], vals2[:,iy],
                plot = result,
                lineStyle = noncomp_style,
                fillcolor = noncomp_color,
                lineTitle = "Non-compatible")
    result.plot(xans, yans, "k-", label="Solution")
    return result

#**************************************** Planar ****************************************
# Read the tabulated point info
labels, vals_comp = readGnuFile("dumps-Sedov-RZ/planar/SPH/nPerh=6.010000/compatibleEnergy/Sedov-planar-RZ.gnu")
labels, vals_noncomp = readGnuFile("dumps-Sedov-RZ/planar/SPH/nPerh=6.010000/nonconservative/Sedov-planar-RZ.gnu")
ix = labels.index("x")

answer = SedovSolution(1, gamma, E0=E0, rho0=rho0)
xmin, xmax = np.min(vals_comp[:,ix]), np.max(vals_comp[:,ix])
xans = np.linspace(xmin, xmax, nsol, endpoint=True)
nu1 = 1.0/(answer.nu + 2.0)
nu2 = 2.0*nu1
tans = (rshock*(answer.alpha*rho0/E0)**nu1)**(1.0/nu2)
xans, vans, uans, rhoans, Pans, Aans, hans = answer.solution(tans, xans)

# Density
rhoPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "rho",
                 xans, rhoans)
plt.xlabel("$z$")
plt.ylabel("$\\rho$")
rhoPlot.axes.legend()
plt.tight_layout()
rhoPlot.figure.savefig("Sedov-RZ-planar-rho.png")

# Pressure
Pplot = plotIt(labels, vals_comp, vals_noncomp,
               "x", "P",
               xans, Pans)
plt.xlabel("$z$")
plt.ylabel("$P$")
plt.tight_layout()
Pplot.figure.savefig("Sedov-RZ-planar-P.png")

# Velocity
velPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "v",
                 xans, vans)
plt.xlabel("$z$")
plt.ylabel("$v$")
plt.tight_layout()
velPlot.figure.savefig("Sedov-RZ-planar-vel.png")

# eps
epsPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "eps",
                 xans, uans)
plt.ylim(-0.1, 5.5)
plt.xlabel("$z$")
plt.ylabel("$u$")
plt.tight_layout()
epsPlot.figure.savefig("Sedov-RZ-planar-eps.png")

#**************************************** Cylindrical ****************************************
# Read the tabulated point info
labels, vals_comp = readGnuFile("dumps-Sedov-RZ/cylindrical/SPH/nPerh=6.010000/compatibleEnergy/Sedov-cylindrical-RZ.gnu")
labels, vals_noncomp = readGnuFile("dumps-Sedov-RZ/cylindrical/SPH/nPerh=6.010000/nonconservative/Sedov-cylindrical-RZ.gnu")
ix = labels.index("x")

answer = SedovSolution(2, gamma, E0=E0, rho0=rho0)
xmin, xmax = np.min(vals_comp[:,ix]), np.max(vals_comp[:,ix])
xans = np.linspace(xmin, xmax, nsol, endpoint=True)
nu1 = 1.0/(answer.nu + 2.0)
nu2 = 2.0*nu1
tans = (rshock*(answer.alpha*rho0/E0)**nu1)**(1.0/nu2)
xans, vans, uans, rhoans, Pans, Aans, hans = answer.solution(tans, xans)

# Density
rhoPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "rho",
                 xans, rhoans)
plt.xlabel("$r$")
plt.ylabel("$\\rho$")
rhoPlot.axes.legend()
plt.tight_layout()
rhoPlot.figure.savefig("Sedov-RZ-cylindrical-rho.png")

# Pressure
Pplot = plotIt(labels, vals_comp, vals_noncomp,
               "x", "P",
               xans, Pans)
plt.xlabel("$r$")
plt.ylabel("$P$")
plt.tight_layout()
Pplot.figure.savefig("Sedov-RZ-cylindrical-P.png")

# Velocity
velPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "v",
                 xans, vans)
plt.xlabel("$r$")
plt.ylabel("$v$")
plt.tight_layout()
velPlot.figure.savefig("Sedov-RZ-cylindrical-vel.png")

# eps
epsPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "eps",
                 xans, uans)
#plt.ylim(-1.0, 60.0)
plt.xlabel("$r$")
plt.ylabel("$u$")
plt.yscale("log")
plt.ylim(1e-2, 100.0)
plt.tight_layout()
epsPlot.figure.savefig("Sedov-RZ-cylindrical-eps.png")

#**************************************** Spherical ****************************************
# Read the tabulated point info
labels, vals_comp = readGnuFile("dumps-Sedov-RZ/spherical/SPH/nPerh=6.010000/compatibleEnergy/Sedov-spherical-RZ.gnu")
labels, vals_noncomp = readGnuFile("dumps-Sedov-RZ/spherical/SPH/nPerh=6.010000/nonconservative/Sedov-spherical-RZ.gnu")
ix = labels.index("x")

answer = SedovSolution(3, gamma, E0=E0, rho0=rho0)
xmin, xmax = np.min(vals_comp[:,ix]), np.max(vals_comp[:,ix])
xans = np.linspace(xmin, xmax, nsol, endpoint=True)
nu1 = 1.0/(answer.nu + 2.0)
nu2 = 2.0*nu1
tans = (rshock*(answer.alpha*rho0/E0)**nu1)**(1.0/nu2)
xans, vans, uans, rhoans, Pans, Aans, hans = answer.solution(tans, xans)

# Density
rhoPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "rho",
                 xans, rhoans)
plt.xlabel("$\\sqrt{r^2 + z^2}$")
plt.ylabel("$\\rho$")
rhoPlot.axes.legend()
plt.tight_layout()
rhoPlot.figure.savefig("Sedov-RZ-spherical-rho.png")

# Pressure
Pplot = plotIt(labels, vals_comp, vals_noncomp,
               "x", "P",
               xans, Pans)
plt.xlabel("$\\sqrt{r^2 + z^2}$")
plt.ylabel("$P$")
plt.tight_layout()
Pplot.figure.savefig("Sedov-RZ-spherical-P.png")

# Velocity
velPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "v",
                 xans, vans)
plt.xlabel("$\\sqrt{r^2 + z^2}$")
plt.ylabel("$v$")
plt.tight_layout()
velPlot.figure.savefig("Sedov-RZ-spherical-vel.png")

# eps
epsPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "eps",
                 xans, uans)
#plt.ylim(-1.0, 60.0)
plt.xlabel("$\\sqrt{r^2 + z^2}$")
plt.ylabel("$u$")
plt.yscale("log")
plt.ylim(1e-2, 500.0)
plt.tight_layout()
epsPlot.figure.savefig("Sedov-RZ-spherical-eps.png")
