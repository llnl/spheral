from Spheral import *
from SpheralTestUtilities import *
#import Gnuplot
import numpy
from SpheralPlotUtilities import *

################################################################################
def plotW(plot, W, xmin=0.0, xmax=2.0, numPnts=200, Hdet=1.0, title='',
          lineTitle=''):
    dx = (xmax - xmin)/(numPnts - 1)
    x = numpy.array(list(range(numPnts)))
    y = numpy.array([0.0]*numPnts)
    i1,i2,i3 = 0,0,0
    x = dx*x + xmin
    for i in range(numPnts):
        y[i] = W(x[i], Hdet)
        if (i>0):
            i1 += abs(y[i]+y[i-1])/2.0 * dx
            i2 += abs(y[i]+y[i-1])/2.0 * dx * x[i]
            i3 += abs(y[i]+y[i-1])/2.0 * dx * x[i]*x[i]
    Pi = 3.14159
    print("{0:3.3f} {1:3.3f} {2:3.3f}".format(i1*2.0,i2*2.0,i3))
    plot('set xrange [%f:%f]' % (xmin, xmax))
    plot.xlabel('r')
    plot.ylabel('W(r)')
    if title:
        plot.title(title)
    data = Gnuplot.Data(x, y, with_='lines', title=lineTitle)
    plot.replot(data)
    return

import sys, string
kernels = list(map(string.lower, sys.argv[1:]))
print(kernels)

numPts = 51
dx = 1.0/(numPts - 1)

################################################################################
numPoints = 100

kernelDict = {'spline': [BSplineKernel1d(),
                         BSplineKernel2d(),
                         BSplineKernel3d()],
              'w4spline': [W4SplineKernel1d(),
                           W4SplineKernel2d(),
                           W4SplineKernel3d()],
              'wendlandc4': [WendlandC4Kernel1d(),
                             WendlandC4Kernel2d(),
                             WendlandC4Kernel3d()],
              'wendlandc6': [WendlandC6Kernel1d(),
                             WendlandC6Kernel2d(),
                             WendlandC6Kernel3d()],
              'expinv': [ExpInvKernel1d(),
                         ExpInvKernel2d(),
                         ExpInvKernel3d()],
              'quartic': [QuarticSplineKernel1d(),
                          QuarticSplineKernel2d(),
                          QuarticSplineKernel3d()],
##              'gauss': [GaussianKernel1d(3),
##                        GaussianKernel2d(3),
##                        GaussianKernel3d(3)],
##              'supergauss': [SuperGaussianKernel1d(),
##                             SuperGaussianKernel2d(),
##                             SuperGaussianKernel3d()],
##              'pigauss': [PiGaussianKernel1d(1.0),
##                          PiGaussianKernel2d(1.0),
##                          PiGaussianKernel3d(1.0)],
##              'sinc': [SincKernel1d(2),
##                       SincKernel2d(2),
##                       SincKernel3d(2)],
##              'poly1': [NSincPolynomialKernel1d(1),
##                        NSincPolynomialKernel2d(1),
##                        NSincPolynomialKernel3d(1)],
##              'poly3': [NSincPolynomialKernel1d(3),
##                        NSincPolynomialKernel2d(3)],
##              'poly5': [NSincPolynomialKernel1d(5),
##                        NSincPolynomialKernel2d(5)],
##              'poly7': [NSincPolynomialKernel1d(7),
##                        NSincPolynomialKernel2d(7)],
              'spline3': [NBSplineKernel1d(3),
                          NBSplineKernel2d(3),
                          NBSplineKernel3d(3)],
              'spline5': [NBSplineKernel1d(5),
                          NBSplineKernel2d(5),
                          NBSplineKernel3d(5)],
              'spline7': [NBSplineKernel1d(7),
                          NBSplineKernel2d(7),
                          NBSplineKernel3d(7)],
              'spline9': [NBSplineKernel1d(9),
                          NBSplineKernel2d(9),
                          NBSplineKernel3d(9)],
##              'spline11': [NBSplineKernel1d(11),
##                           NBSplineKernel2d(11),
##                           NBSplineKernel3d(11)],
              }

titleDict = {'spline':  'B Spline Kernel',
             'h': 'H kernel',
             'h10': 'H kernel (extent = 10)',
             'quartic': 'Quartic Spline Kernel',
             'w4spline': 'W4 Spline Kernel',
             'gauss': 'Gaussian Kernel',
             'supergauss': 'SuperGaussian Kernel',
             'pigauss': 'Pi Gaussian Kernel',
             'sinc': 'Sinc Kernel',
             'poly1': 'Linear Polynomial Sinc approx Kernel',
             'poly3': 'Cubic Polynomial Sinc approx Kernel',
             'poly5': 'Quintic Polynomial Sinc approx Kernel',
             'poly7': 'Septic Polynomial Sinc approx Kernel',
             'spline3': '3rd order b spline Kernel',
             'spline5': '5th order b spline Kernel',
             'spline7': '7th order b spline Kernel',
             'spline9': '9th order b spline Kernel',
             'spline11': '11th order b spline Kernel',
             'wendlandc4': 'Wendland C4 Kernel',
             'wendlandc6': 'Wendland C6 Kernel',
             'expinv' : 'Exponential inverse Kernel',
             }

plots = []
for kernel in kernels:
    title(titleDict[kernel])
    for i in range(6):
        plots.append(generateNewGnuPlot())
    for W in kernelDict[kernel]:
        output("W")
        output("W.volumeNormalization")
        output("W.kernelExtent")
        plotW(plots[-6], W.kernelValue, 0.0, W.kernelExtent,
              title = titleDict[kernel],
              lineTitle = str(W))
        plotW(plots[-5], W.gradValue, 0.0, W.kernelExtent,
              title = titleDict[kernel] + ' Gradient',
              lineTitle = str(W))
        plotW(plots[-4], W.grad2Value, 0.0, W.kernelExtent,
              title = titleDict[kernel] + ' second derivative',
              lineTitle = str(W))

        # Build a tabular version of the kernel
        WT = eval('TableKernel' + str(W).split()[0][-2:] + '(W, numPoints)')
        output("WT")
        output("WT.volumeNormalization")
        output("WT.kernelExtent")
        plotW(plots[-3], WT.kernelValue, 0.0, W.kernelExtent,
              title='Table ' + titleDict[kernel],
              lineTitle = str(W))
        plotW(plots[-2], WT.gradValue, 0.0, W.kernelExtent,
              title = 'Table ' + titleDict[kernel] + ' Gradient',
              lineTitle = str(W))
        plotW(plots[-1], WT.grad2Value, 0.0, W.kernelExtent,
              title = 'Table ' + titleDict[kernel] + ' second derivative',
              lineTitle = str(W))

