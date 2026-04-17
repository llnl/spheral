#-------------------------------------------------------------------------------
# Analysis methods to plot the Noh RZ results
#-------------------------------------------------------------------------------
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import numpy as np
from NohAnalyticSolution import NohSolution

plt.rc("font", size=18)

# Parameters
gamma = 5.0/3.0
rho0 = 1.0
nsol = 1000
tans = 0.6
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
labels, vals_comp = readGnuFile("dumps-rz-Noh/planar/SPH/nPerh=6.010000/compatibleEnergy=True/Cullen=False/Noh-RZ.gnu")
labels, vals_noncomp = readGnuFile("dumps-rz-Noh/planar/SPH/nPerh=6.010000/compatibleEnergy=False/Cullen=False/Noh-RZ.gnu")
ix = labels.index("x")

xmin, xmax = np.min(vals_comp[:,ix]), np.max(vals_comp[:,ix])
xans = np.linspace(xmin, xmax, nsol, endpoint=True)
answer = NohSolution(1, r=xans, gamma=gamma, rho0=rho0)
xans, vans, uans, rhoans, Pans, hans = answer.solution(tans, xans)

# Density
rhoPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "rho",
                 xans, rhoans)
plt.ylim(0.5, 4.5)
plt.xlabel("$z$")
plt.ylabel("$\\rho$")
plt.tight_layout()
rhoPlot.figure.savefig("Noh-RZ-planar-rho.png")

# Pressure
Pplot = plotIt(labels, vals_comp, vals_noncomp,
               "x", "P",
               xans, Pans)
plt.xlabel("$z$")
plt.ylabel("$P$")
plt.tight_layout()
Pplot.figure.savefig("Noh-RZ-planar-P.png")

# Velocity
velPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "v",
                 xans, vans)
plt.xlabel("$z$")
plt.ylabel("$v$")
plt.tight_layout()
velPlot.figure.savefig("Noh-RZ-planar-vel.png")

# eps
epsPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "eps",
                 xans, uans)
plt.ylim(-0.05, 0.6)
plt.xlabel("$z$")
plt.ylabel("$u$")
epsPlot.axes.legend()
plt.tight_layout()
epsPlot.figure.savefig("Noh-RZ-planar-eps.png")

#**************************************** Cylindrical ****************************************
# Read the tabulated point info
labels, vals_comp = readGnuFile("dumps-rz-Noh/cylindrical/SPH/nPerh=6.010000/compatibleEnergy=True/Cullen=False/Noh-RZ.gnu")
labels, vals_noncomp = readGnuFile("dumps-rz-Noh/cylindrical/SPH/nPerh=6.010000/compatibleEnergy=False/Cullen=False/Noh-RZ.gnu")
ix = labels.index("x")

xmin, xmax = np.min(vals_comp[:,ix]), np.max(vals_comp[:,ix])
xans = np.linspace(xmin, xmax, nsol, endpoint=True)
answer = NohSolution(2, r=xans, gamma=gamma, rho0=rho0)
xans, vans, uans, rhoans, Pans, hans = answer.solution(tans, xans)

# Density
rhoPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "rho",
                 xans, rhoans)
plt.xlabel("$r$")
plt.ylabel("$\\rho$")
plt.tight_layout()
rhoPlot.figure.savefig("Noh-RZ-cylindrical-rho.png")

# Pressure
Pplot = plotIt(labels, vals_comp, vals_noncomp,
               "x", "P",
               xans, Pans)
plt.xlabel("$r$")
plt.ylabel("$P$")
plt.tight_layout()
Pplot.figure.savefig("Noh-RZ-cylindrical-P.png")

# Velocity
velPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "v",
                 xans, vans)
plt.xlabel("$r$")
plt.ylabel("$v$")
plt.tight_layout()
velPlot.figure.savefig("Noh-RZ-cylindrical-vel.png")

# eps
epsPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "eps",
                 xans, uans)
plt.ylim(-0.05, 1.0)
plt.xlabel("$r$")
plt.ylabel("$u$")
epsPlot.axes.legend()
plt.tight_layout()
epsPlot.figure.savefig("Noh-RZ-cylindrical-eps.png")

#**************************************** Spherical ****************************************
# Read the tabulated point info
labels, vals_comp = readGnuFile("dumps-rz-Noh/spherical/SPH/nPerh=6.010000/compatibleEnergy=True/Cullen=False/Noh-RZ.gnu")
labels, vals_noncomp = readGnuFile("dumps-rz-Noh/spherical/SPH/nPerh=6.010000/compatibleEnergy=False/Cullen=False/Noh-RZ.gnu")
ix = labels.index("x")

xmin, xmax = np.min(vals_comp[:,ix]), np.max(vals_comp[:,ix])
xans = np.linspace(xmin, xmax, nsol, endpoint=True)
answer = NohSolution(3, r=xans, gamma=gamma, rho0=rho0)
xans, vans, uans, rhoans, Pans, hans = answer.solution(tans, xans)

# Density
rhoPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "rho",
                 xans, rhoans)
plt.xlabel("$\\sqrt{r^2 + z^2}$")
plt.ylabel("$\\rho$")
plt.tight_layout()
rhoPlot.figure.savefig("Noh-RZ-spherical-rho.png")

# Pressure
Pplot = plotIt(labels, vals_comp, vals_noncomp,
               "x", "P",
               xans, Pans)
plt.xlabel("$\\sqrt{r^2 + z^2}$")
plt.ylabel("$P$")
plt.tight_layout()
Pplot.figure.savefig("Noh-RZ-spherical-P.png")

# Velocity
velPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "v",
                 xans, vans)
plt.xlabel("$\\sqrt{r^2 + z^2}$")
plt.ylabel("$v$")
plt.tight_layout()
velPlot.figure.savefig("Noh-RZ-spherical-vel.png")

# eps
epsPlot = plotIt(labels, vals_comp, vals_noncomp,
                 "x", "eps",
                 xans, uans)
plt.ylim(-0.05, 1.0)
plt.xlabel("$\\sqrt{r^2 + z^2}$")
plt.ylabel("$u$")
epsPlot.axes.legend()
plt.tight_layout()
epsPlot.figure.savefig("Noh-RZ-spherical-eps.png")
