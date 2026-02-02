import numpy as np
import sys

from Spheral import *
from SpheralTestUtilities import *
from SpheralMatplotlib import *

#-------------------------------------------------------------------------------
# What kernels should we plot
#-------------------------------------------------------------------------------
kernels = sys.argv[1:]
output("kernels")

#-------------------------------------------------------------------------------
# SPH zeroth moment algorithm
#-------------------------------------------------------------------------------
def sumKernelValues1d(WT, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = sum([abs(WT.kernelValueSPH(abs(etax))) for etax in np.arange(-etamax, etamax, deta)])
    return result

def sumKernelValues2d(WT, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = 0.0
    for etay in np.arange(-etamax, etamax, deta):
        for etax in np.arange(-etamax, etamax, deta):
            eta = sqrt(etax*etax + etay*etay)
            result += WT.kernelValueSPH(eta)
    return sqrt(result)

def sumKernelValues3d(WT, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = 0.0
    for etaz in np.arange(-etamax, etamax, deta):
        for etay in np.arange(-etamax, etamax, deta):
            for etax in np.arange(-etamax, etamax, deta):
                eta = sqrt(etax*etax + etay*etay + etaz*etaz)
                result += WT.kernelValueSPH(eta)
    return (result)**(1.0/3.0)

#-------------------------------------------------------------------------------
# ASPH second moment algorithm
#-------------------------------------------------------------------------------
def safeInv(x, fuzz=1e-30):
    return x/(x*x + fuzz)

def sumKernelValuesASPH1d(WT, targetNperh, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = sum([WT.kernelValueASPH(abs(etax), targetNperh)*etax*etax for etax in np.arange(-etamax, etamax, deta)])
    return result

def sumKernelValuesASPH2d(WT, targetNperh, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    result = SymTensor2d()
    for etay in np.arange(-etamax, etamax, deta):
        for etax in np.arange(-etamax, etamax, deta):
            eta = Vector2d(etax, etay)
            Wi = WT.kernelValueASPH(eta.magnitude(), targetNperh)
            result += Wi * eta.selfdyad()
    return sqrt(0.5*result.eigenValues().sumElements())

def sumKernelValuesASPH3d(WT, targetNperh, nperh):
    deta = 1.0/nperh
    etamax = WT.kernelExtent
    Wsum = SymTensor3d()
    result = SymTensor3d()
    for etaz in np.arange(-etamax, etamax, deta):
        for etay in np.arange(-etamax, etamax, deta):
            for etax in np.arange(-etamax, etamax, deta):
                eta = Vector3d(etax, etay, etaz)
                result += WT.kernelValueASPH(eta.magnitude(), targetNperh) * eta.selfdyad()
    return (result.eigenValues().sumElements()/3.0)**(1.0/3.0)

# def sumKernelValuesSlice2d(WT, nhat, nperh, detax, detay):
#     etamax = WT.kernelExtent
#     result = SymTensor2d()
#     for etay in np.arange(-etamax, etamax, detay):
#         for etax in np.arange(-etamax, etamax, detax):
#             eta = Vector2d(etax, etay)
#             result += WT.kernelValueASPH(eta.magnitude(), nperh) * eta.selfdyad()
#     return sqrt((result*nhat).magnitude())

# def sumKernelValuesSlice3d(WT, nhat, nperh, detax, detay, detaz):
#     etamax = WT.kernelExtent
#     result = SymTensor3d()
#     for etaz in np.arange(-etamax, etamax, detaz):
#         for etay in np.arange(-etamax, etamax, detay):
#             for etax in np.arange(-etamax, etamax, detax):
#                 eta = Vector3d(etax, etay, etaz)
#                 result += WT.kernelValueASPH(eta.magnitude(), nperh) * eta.selfdyad()
#     return ((result*nhat).magnitude())**(1.0/3.0)

#-------------------------------------------------------------------------------
# Here we go...
#-------------------------------------------------------------------------------
for Wstr in kernels:
    title(Wstr)

    nDim = 0
    if str(Wstr).split()[0][-2:] == "1d":
        nDim = 1
    elif str(Wstr).split()[0][-2:] == "2d":
        nDim = 2
    elif str(Wstr).split()[0][-2:] == "3d":
        nDim = 3
    assert nDim > 0

    # Plot the kernel basics
    WT = eval(f"TableKernel{nDim}d({Wstr}(), 400)")
    #plotTableKernel(WT, nPerh=4.01)

    targetNperh = 4.01
    asph = eval(f"ASPHSmoothingScale{nDim}d(WT, {targetNperh}, 400)")

    # Now how well do we recover nPerh based on kernel sums?
    etamax = WT.kernelExtent
    nperh0 = np.arange(1.0/etamax, 10.0, 0.1)
    nperhSPH = []
    nperhASPH = []
    WsumSPH = []
    WsumASPH = []
    for nperh in nperh0:
        Wsumi = eval(f"sumKernelValues{nDim}d(WT, {nperh})")
        WsumASPHi = eval(f"sumKernelValuesASPH{nDim}d(WT, {targetNperh}, {nperh})")
        WsumSPH.append(Wsumi)
        WsumASPH.append(WsumASPHi)
        nperhSPH.append(WT.equivalentNodesPerSmoothingScale(Wsumi))
        nperhASPH.append(asph.equivalentNodesPerSmoothingScale(WsumASPHi))
    WsumSPH = np.array(WsumSPH)
    WsumASPH = np.array(WsumASPH)
    nperhSPH = np.array(nperhSPH)
    nperhASPH = np.array(nperhASPH)

    # Helper function for plotting
    def plotIt(x, y, style,
               label = None,
               xlabel = None,
               ylabel = None,
               title = None,
               plot = None):
        if plot is None:
            plot = newFigure()
        plot.plot(x, y, style, label=label)
        if title:
            plot.set_title(title)
        if xlabel:
            plot.set_xlabel(xlabel)
        if ylabel:
            plot.set_ylabel(ylabel)
        plot.legend()
        return plot

    # SPH fit for nperh(Wsum)
    plot = plotIt(WsumSPH, nperh0, "r-*", label="Actual",
                  title = f"{Wstr} n per h as a function of $\sum W$ : SPH algorithm",
                  xlabel = r"$\sum W$",
                  ylabel = "n per h")
    plotIt(WsumSPH, nperhSPH, "k-", label="Fit", plot=plot)

    # ASPH fit for nperh(Wsum)
    plot = plotIt(WsumASPH, nperh0, "r-*", label="Actual",
                  title = f"{Wstr} n per h as a function of $\lambda(\psi)$ : ASPH algorithm",
                  xlabel = r"$\lambda(\psi)$",
                  ylabel = "n per h")
    plotIt(WsumASPH, nperhASPH, "k-", label="Fit", plot=plot)

    # # SPH nperh
    # plot = plotIt(nperh0, nperhSPH, "b*-", label="nperh lookup",
    #               title = f"{Wstr} n per h lookup test : SPH algorithm",
    #               xlabel = "nperh actual",
    #               ylabel = "nperh estimated")

    # # SPH nperh error
    # plot = plotIt(nperh0, (nperhSPH - nperh0)/nperh0, "r*-",
    #               title = f"{Wstr} n per h lookup test error : SPH algorithm",
    #               xlabel = "nperh actual",
    #               ylabel = "Error")

    # plot = plotIt(nperh0, nperhASPH, "b*-",
    #               title = f"{Wstr} n per h lookup test : ASPH algorithm",
    #               xlabel = "nperh actual",
    #               ylabel = "nperh estimated")

    # plot = plotIt(nperh0, (nperhASPH - nperh0)/nperh0, "r*-",
    #               title = f"{Wstr} n per h lookup test error : ASPH algorithm",
    #               xlabel = "nperh actual",
    #               ylabel = "Error")

    # # Test ASPH with different aspect ratios
    # if nDim == 2:
    #     aspect = np.arange(0.1, 1.0, 0.05)
    #     X, Y = np.meshgrid(nperh0, aspect)
    #     WsumASPHx = np.ndarray(X.shape)
    #     WsumASPHy = np.ndarray(X.shape)
    #     nperhASPHx = np.ndarray(X.shape)
    #     nperhASPHy = np.ndarray(X.shape)
    #     nperhASPHx_err = np.ndarray(X.shape)
    #     nperhASPHy_err = np.ndarray(X.shape)
    #     for iy in range(X.shape[0]):
    #         for ix in range(X.shape[1]):
    #             nPerhi = X[iy,ix]
    #             aspecti = Y[iy,ix]
    #             WsumASPHx[iy,ix] = sumKernelValuesSlice2d(WT, Vector2d(1,0), nPerhi, 1.0/nPerhi, aspecti/nPerhi)
    #             WsumASPHy[iy,ix] = sumKernelValuesSlice2d(WT, Vector2d(0,1), nPerhi, 1.0/nPerhi, aspecti/nPerhi)
    #             nperhASPHx[iy,ix] = WT.equivalentNodesPerSmoothingScaleASPH(WsumASPHx[iy,ix])
    #             nperhASPHy[iy,ix] = WT.equivalentNodesPerSmoothingScaleASPH(WsumASPHy[iy,ix])
    #             nperhASPHx_err[iy,ix] = (nperhASPHx[iy,ix] - nPerhi)/nPerhi
    #             nperhASPHy_err[iy,ix] = (nperhASPHy[iy,ix] - nPerhi/aspecti)/(nPerhi/aspecti)

    #     plotSurface(X, Y, WsumASPHx,
    #                 title = f"{Wstr} ASPH Wsum $X$",
    #                 xlabel = "n per h",
    #                 ylabel = "aspect ratio")
    #     plotSurface(X, Y, WsumASPHy,
    #                 title = f"{Wstr} ASPH Wsum $Y$",
    #                 xlabel = "n per h",
    #                 ylabel = "aspect ratio")
    #     plotSurface(X, Y, nperhASPHx,
    #                 title = f"{Wstr} ASPH n per h $X$",
    #                 xlabel = "n per h",
    #                 ylabel = "aspect ratio")
    #     plotSurface(X, Y, nperhASPHy,
    #                 title = f"{Wstr} ASPH n per h $Y$",
    #                 xlabel = "n per h",
    #                 ylabel = "aspect ratio")
    #     plotSurface(X, Y, nperhASPHx_err,
    #                 title = f"{Wstr} ASPH n per h $X$ error",
    #                 xlabel = "n per h",
    #                 ylabel = "aspect ratio")
    #     plotSurface(X, Y, nperhASPHy_err,
    #                 title = f"{Wstr} ASPH n per h $Y$ error",
    #                 xlabel = "n per h",
    #                 ylabel = "aspect ratio")
