from matplotlib.pyplot import cm as pltcm
from matplotlib import patches
#from matplotlib.collections import PatchCollections
import numpy as np
import mpi
from Spheral import *
from math import *
import os
from SpheralTestUtilities import multiSort

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Dummy type to handle parallel matplotlib
#-------------------------------------------------------------------------------
class NullFigure:
    def __init__(self):
        pass
    def __call__(self, *arghs, **keyw):
        pass
    def __getattr__(self, name, *args, **keyw):
        def method(*args, **keyw):
            pass
        if args or keyw:
            return method
        else:
            return NullFigure()
    def __setattr__(self, name, val):
        pass
    def savefig(self, *arghs, **keyw):
        pass

#-------------------------------------------------------------------------------
# Parallel safe pyplot
#-------------------------------------------------------------------------------
if mpi.rank == 0:
    import matplotlib.pyplot as plt
else:
    plt = NullFigure()

# Turn on TeX
plt.rc("text", usetex = True)

#-------------------------------------------------------------------------------
# Parallel safe new plot.
#-------------------------------------------------------------------------------
def newFigure():
    if mpi.rank == 0:
        return plt.figure().add_subplot(111)
    else:
        return NullFigure()

#-------------------------------------------------------------------------------
# Take one of our plotted objects and save it to a file.
#-------------------------------------------------------------------------------
def savefig(plot,
            fname,
            transparent = False):
    plot.figure.savefig(fname, transparent=transparent, bbox_inches="tight")

#-------------------------------------------------------------------------------
# Calculate the radial velocity component, given a FieldList of positions
# and a FieldList of velocities.
#-------------------------------------------------------------------------------
def radialVelocityFieldList(positions,
                            velocities):

    dim = type(positions).__name__[-2:]
    radialVelocity = None
    fieldConstructor = None
    if dim == "1d":
        radialVelocity = ScalarFieldList1d()
        fieldConstructor = ScalarField1d
    elif dim == "2d":
        radialVelocity = ScalarFieldList2d()
        fieldConstructor = ScalarField2d
    elif dim == "3d":
        radialVelocity = ScalarFieldList3d()
        fieldConstructor = ScalarField3d
    radialVelocity.copyFields()
    for field in positions:
        radialVelocity.appendField(fieldConstructor("radial velocity", field.nodeList()))

    assert positions.numFields == velocities.numFields == radialVelocity.numFields

    for fieldID in range(positions.numFields):
        rfield = positions[fieldID]
        vfield = velocities[fieldID]
        vrfield = radialVelocity[fieldID]
        assert rfield.numElements == vfield.numElements == vrfield.numElements

        for nodeID in range(rfield.numElements):
            r = rfield[nodeID]
            v = vfield[nodeID]
            runit = r.unitVector()
            vrfield[nodeID] = v.dot(runit)

    return radialVelocity

#-------------------------------------------------------------------------------
# Calculate the azimuthal velocity component, given a FieldList of positions
# and a FieldList of velocities.
#-------------------------------------------------------------------------------
def azimuthalVelocityFieldList(positions,
                               velocities):

    dim = type(positions).__name__[-2:]
    azimuthalVelocity = None
    fieldConstructor = None
    if dim == "1d":
        azimuthalVelocity = ScalarFieldList1d()
        fieldConstructor = ScalarField1d
    elif dim == "2d":
        azimuthalVelocity = ScalarFieldList2d()
        fieldConstructor = ScalarField2d
    elif dim == "3d":
        azimuthalVelocity = ScalarFieldList3d()
        fieldConstructor = ScalarField3d
    azimuthalVelocity.copyFields()
    for field in positions:
        azimuthalVelocity.appendField(fieldConstructor("azimuthal velocity", field.nodeList()))

    assert positions.numFields == velocities.numFields == azimuthalVelocity.numFields

    for fieldID in range(positions.numFields):
        rfield = positions[fieldID]
        vfield = velocities[fieldID]
        vafield = azimuthalVelocity[fieldID]
        assert rfield.numElements == vfield.numElements == vafield.numElements

        for nodeID in range(rfield.numElements):
            r = rfield[nodeID]
            v = vfield[nodeID]
            raz = r.unitVector()
            x = raz.x
            y = raz.y
            raz.x = -y
            raz.y = x
            vafield[nodeID] = v.dot(raz)

    return azimuthalVelocity

#-------------------------------------------------------------------------------
# Helper method to determine the angular momentum per node.
#-------------------------------------------------------------------------------
def angularMomentum(mass, position, velocity):
    assert mass.numFields == position.numFields == velocity.numFields

    result = []
    for massField, positionField, velocityField in zip(mass,
                                                       position,
                                                       velocity):
        assert (massField.nodeList().numInternalNodes ==
                positionField.nodeList().numInternalNodes ==
                velocityField.nodeList().numInternalNodes)
        for j in range(massField.nodeList().numInternalNodes):
            result.append((positionField[j].cross(velocityField[j]))*massField[j])

    return result

#-------------------------------------------------------------------------------
# Plot a FieldList
#-------------------------------------------------------------------------------
def plotFieldList(fieldList,
                  xFunction = "%s.x",
                  yFunction = "%s",
                  plotGhosts = False,
                  colorNodeLists = False,
                  plot = None,
                  xRange = [None, None],
                  yRange = [None, None],
                  plotStyle = "ro",
                  markerSize = 4,
                  kwords = {},
                  winTitle = None,
                  lineTitle = "",
                  xlabel = None,
                  ylabel = None,
                  filterFunc = None,
                  semilogy = False):

    # Do we need to make a new window?
    if plot is None:
        plot = newFigure()

    # How about a filtering function?
    if filterFunc is None:
        filterFunc = lambda x: True

    # Gather the fieldList info across all processors to process 0.
    globalNumNodes = []
    globalX = []
    globalY = []
    for field in fieldList:
        if plotGhosts:
            xvals = field.nodeList().positions().allValues()
            yvals = field.allValues()
        else:
            xvals = field.nodeList().positions().internalValues()
            yvals = field.internalValues()
        localX = []
        localY = []
        for x, y in zip(xvals, yvals):
            if filterFunc(x):
                localX.append(eval(xFunction % "x"))
                localY.append(eval(yFunction % "y"))
        n = len(localX)
        globalNumNodes.append(mpi.allreduce(n, mpi.SUM))
        globalX += mpi.allreduce(localX, mpi.SUM)
        globalY += mpi.allreduce(localY, mpi.SUM)
            
    if mpi.rank == 0:
        # Find the total number of nodes.
        totalNumNodes = sum(globalNumNodes)
        assert(len(globalNumNodes) == fieldList.numFields)
        assert(len(globalX) == totalNumNodes)
        assert(len(globalY) == totalNumNodes)

        # Set the labels.
        if winTitle: plt.title(winTitle)
        if xlabel: plt.xlabel(xlabel)
        if ylabel: plt.ylabel(ylabel)

        # Finally, loop over the fields and do the deed.
        assert len(globalX) == len(globalY)
        if colorNodeLists:
            cumulativeNumNodes = 0
            for fieldID in range(len(globalNumNodes)):
                n = globalNumNodes[fieldID]
                if n:
                    if semilogy:
                        plot.semilogy(globalX[cumulativeNumNodes:cumulativeNumNodes + n],
                                      globalY[cumulativeNumNodes:cumulativeNumNodes + n],
                                      plotStyle, ms=markerSize, label = "%s: %s" % (lineTitle, fieldList[i].nodeList().name), **kwords)
                    else:
                        plot.plot(globalX[cumulativeNumNodes:cumulativeNumNodes + n],
                                  globalY[cumulativeNumNodes:cumulativeNumNodes + n],
                                  plotStyle, ms=markerSize, label = "%s: %s" % (lineTitle, fieldList[i].nodeList().name), **kwords)
                    cumulativeNumNodes += n
        else:
            if semilogy:
                plot.semilogy(globalX, globalY, plotStyle, ms=markerSize, label = lineTitle, **kwords)
            else:
                plot.plot(globalX, globalY, plotStyle, ms=markerSize, label = lineTitle, **kwords)
        plot.axes.legend()

        # Set the ranges.
        xmin, xmax = plt.xlim()
        ymin, ymax = plt.ylim()
        if xRange[0]: xmin = xRange[0]
        if xRange[1]: xmax = xRange[1]
        if yRange[0]: ymin = yRange[0]
        if yRange[1]: ymax = yRange[1]
        plt.xlim(xmin, xmax)
        plt.ylim(ymin, ymax)

    # That's it
    mpi.barrier()
    return plot

#-------------------------------------------------------------------------------
# Plot a Field by promoting it to a FieldList
#-------------------------------------------------------------------------------
def plotField(field,
              xFunction = "%s.x",
              yFunction = "%s",
              plotGhosts = False,
              colorNodeLists = False,
              plot = None,
              xRange = [None, None],
              yRange = [None, None],
              plotStyle = "ro",
              markerSize = 4,
              kwords = {},
              winTitle = None,
              lineTitle = "",
              xlabel = None,
              ylabel = None,
              filterFunc = None,
              semilogy = False):
    flt = eval(type(field).__name__.replace("Field", "FieldList"))
    fl = flt()
    fl.appendField(field)
    return plotFieldList(fl,
                         xFunction = xFunction,
                         yFunction = yFunction,
                         plotGhosts = plotGhosts,
                         colorNodeLists = colorNodeLists,
                         plot = plot,
                         xRange = xRange,
                         yRange = yRange,
                         plotStyle = plotStyle,
                         markerSize = markerSize,
                         kwords = kwords,
                         winTitle = winTitle,
                         lineTitle = lineTitle,
                         xlabel = xlabel,
                         ylabel = ylabel,
                         filterFunc = filterFunc,
                         semilogy = semilogy)

#-------------------------------------------------------------------------------
# Plot the mass density, velocity, pressure, and smoothing scale for the fluid
# node lists in the given data base.  Implicitly assuming 1-D.
#-------------------------------------------------------------------------------
def plotState(thingus,
              plotGhosts = False,
              colorNodeLists = False,
              plotStyle = "ro",
              markerSize = 4,
              xFunction = "%s.x",
              vecyFunction = "%s.x",
              tenyFunction = "%s.xx ** -1",
              lineTitle = "Simulation",
              filterFunc = None):

    dim = type(thingus).__name__[-2:]
    if isinstance(thingus, eval("State%s" % dim)):
        rho = thingus.scalarFields(HydroFieldNames.massDensity)
        vel = thingus.vectorFields(HydroFieldNames.velocity)
        eps = thingus.scalarFields(HydroFieldNames.specificThermalEnergy)
        P = thingus.scalarFields(HydroFieldNames.pressure)
        H = thingus.symTensorFields(HydroFieldNames.H)

    else:
        assert isinstance(thingus, eval("DataBase%s" % dim))
        rho = thingus.fluidMassDensity
        vel = thingus.fluidVelocity
        eps = thingus.fluidSpecificThermalEnergy
        P = thingus.newFluidScalarFieldList(0.0, "pressure")
        thingus.fluidPressure(P)
        H = thingus.fluidHfield

    rhoPlot = plotFieldList(rho,
                            xFunction = xFunction,
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = plotStyle,
                            markerSize = markerSize,
                            winTitle = "Mass Density",
                            lineTitle = lineTitle,
                            xlabel="x",
                            filterFunc = filterFunc)

    velPlot = plotFieldList(vel,
                            xFunction = xFunction,
                            yFunction = vecyFunction,
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = plotStyle,
                            markerSize = markerSize,
                            winTitle = "Velocity",
                            lineTitle = lineTitle,
                            xlabel="x",
                            filterFunc = filterFunc)

    epsPlot = plotFieldList(eps,
                            xFunction = xFunction,
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = plotStyle,
                            markerSize = markerSize,
                            winTitle = "Specific Thermal Energy",
                            lineTitle = lineTitle,
                            xlabel="x",
                            filterFunc = filterFunc)

    PPlot = plotFieldList(P,
                          xFunction = xFunction,
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          plotStyle = plotStyle,
                          markerSize = markerSize,
                          winTitle = "Pressure",
                          lineTitle = lineTitle,
                          xlabel="x",
                          filterFunc = filterFunc)

    HPlot = plotFieldList(H,
                          xFunction = xFunction,
                          yFunction = tenyFunction,
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          plotStyle = plotStyle,
                          markerSize = markerSize,
                          winTitle = "Smoothing scale",
                          lineTitle = lineTitle,
                          xlabel="x",
                          filterFunc = filterFunc)

    return rhoPlot, velPlot, epsPlot, PPlot, HPlot

#-------------------------------------------------------------------------------
# Plot the state vs. radius
#-------------------------------------------------------------------------------
def plotRadialState(dataBase,
                    plotGhosts = False,
                    colorNodeLists = False,
                    lineTitle = "Simulation",
                    filterFunc = None):

    rhoPlot = plotFieldList(dataBase.fluidMassDensity,
                            xFunction = "%s.magnitude()",
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = "ro",
                            winTitle = "Mass density",
                            lineTitle = lineTitle,
                            xlabel = "r",
                            filterFunc = filterFunc)

    radialVelocity = radialVelocityFieldList(dataBase.fluidPosition,
                                             dataBase.fluidVelocity)
    velPlot = plotFieldList(radialVelocity,
                            xFunction = "%s.magnitude()",
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = "ro",
                            winTitle = " Radial Velocity",
                            lineTitle = lineTitle,
                            xlabel = "r",
                            filterFunc = filterFunc)

    epsPlot = plotFieldList(dataBase.fluidSpecificThermalEnergy,
                            xFunction = "%s.magnitude()",
                            plotGhosts = plotGhosts,
                            colorNodeLists = colorNodeLists,
                            plotStyle = "ro",
                            winTitle = "Specific Thermal Energy",
                            lineTitle = lineTitle,
                            xlabel = "r",
                            filterFunc = filterFunc)

    fluidPressure = dataBase.newFluidScalarFieldList(0.0, "pressure")
    dataBase.fluidPressure(fluidPressure)
    PPlot = plotFieldList(fluidPressure,
                          xFunction = "%s.magnitude()",
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          plotStyle = "ro",
                          winTitle = "Pressure",
                          lineTitle = lineTitle,
                          xlabel = "r",
                          filterFunc = filterFunc)

    HPlot = plotFieldList(dataBase.fluidHfield,
                          xFunction = "%s.magnitude()",
                          yFunction = "%s.xx**-1",
                          plotGhosts = plotGhosts,
                          colorNodeLists = colorNodeLists,
                          plotStyle = "ro",
                          winTitle = "Smoothing scale",
                          lineTitle = lineTitle,
                          xlabel = "r",
                          filterFunc = filterFunc)

    return rhoPlot, velPlot, epsPlot, PPlot, HPlot

#-------------------------------------------------------------------------------
# Overplot the answer on results from plotState.
#-------------------------------------------------------------------------------
def plotAnswer(answerObject, time,
               rhoPlot = None,
               velPlot = None,
               epsPlot = None,
               PPlot = None,
               APlot = None,
               HPlot = None,
               x = None,
               plotStyle = "k-"):

    try:
        x, v, u, rho, P, h = answerObject.solution(time, x)
    except:
        try:
            x, v, u, rho, P, A, h = answerObject.solution(time, x)
        except:
            x, v, u, rho, P = answerObject.solution(time, x)

    if rhoPlot is not None:
        rhoPlot.plot(x, rho, plotStyle, label="Solution")
        rhoPlot.axes.legend()

    if velPlot is not None:
        velPlot.plot(x, v, plotStyle, label="Solution")
        velPlot.axes.legend()

    if epsPlot is not None:
        epsPlot.plot(x, u, plotStyle, label="Solution")
        epsPlot.axes.legend()

    if PPlot is not None:
        PPlot.plot(x, P, plotStyle, label="Solution")
        PPlot.axes.legend()

    if APlot is not None:
        APlot.plot(x, A, plotStyle, label="Solution")
        APlot.axes.legend()

    if HPlot is not None:
        HPlot.plot(x, h, plotStyle, label="Solution")
        HPlot.axes.legend()

    return

#-------------------------------------------------------------------------------
# Plot the node positions
#-------------------------------------------------------------------------------
def plotNodePositions2d(thingy,
                        xFunction = "%s.x",
                        yFunction = "%s.y",
                        plotGhosts = False,
                        colorNodeLists = True,
                        colorDomains = False,
                        title = "",
                        plotStyle = "ro",
                        markerSize = 4):

    assert colorNodeLists + colorDomains <= 1

    if isinstance(thingy, DataBase2d):
        nodeLists = thingy.nodeLists
    else:
        nodeLists = thingy

    # Gather the node positions across all domains.
    # Loop over all the NodeLists.
    xNodes = []
    yNodes = []
    for nodeList in nodeLists:
        if plotGhosts:
            pos = nodeList.positions().allValues()
        else:
            pos = nodeList.positions().internalValues()
        xNodes.append([eval(xFunction % "x") for x in pos])
        yNodes.append([eval(yFunction % "x") for x in pos])
    assert len(xNodes) == len(nodeLists)
    assert len(xNodes) == len(yNodes)
    
    globalXNodes = mpi.gather(xNodes)
    globalYNodes = mpi.gather(yNodes)

    plot = newFigure()

    if mpi.rank == 0:
        assert len(globalXNodes) == mpi.procs
        assert len(globalYNodes) == mpi.procs

        xlist, ylist = [], []
        if colorDomains:
            for xDomain, yDomain in zip(globalXNodes, globalYNodes):
                assert len(xDomain) == len(nodeLists)
                assert len(yDomain) == len(nodeLists)
                xlist.append([])
                ylist.append([])
                for xx in xDomain:
                    xlist[-1].extend(xx)
                for yy in yDomain:
                    ylist[-1].extend(yy)
            assert len(xlist) == mpi.procs
            assert len(ylist) == mpi.procs

        elif colorNodeLists:
            for i in range(len(nodeLists)):
                xlist.append([])
                ylist.append([])
            for xDomain, yDomain in zip(globalXNodes, globalYNodes):
                assert len(xDomain) == len(nodeLists)
                assert len(yDomain) == len(nodeLists)
                for i in range(len(nodeLists)):
                    xlist[i].extend(xDomain[i])
                    ylist[i].extend(yDomain[i])
            assert len(xlist) == len(nodeLists)
            assert len(ylist) == len(nodeLists)

        else:
            xlist, ylist = [[]], [[]]
            for xDomain, yDomain in zip(globalXNodes, globalYNodes):
                print(len(xDomain), len(nodeLists))
                assert len(xDomain) == len(nodeLists)
                assert len(yDomain) == len(nodeLists)
                for i in range(len(nodeLists)):
                    xlist[0].extend(xDomain[i])
                    ylist[0].extend(yDomain[i])

        plt.title(title)
        color = iter(pltcm.rainbow(np.linspace(0,1,len(xlist))))
        for x, y in zip(xlist, ylist):
            c = next(color)
            plot.plot(x, y, "o", color=c, ms=markerSize)
        plot.axes.set_aspect("equal", "datalim")

    return plot

#-------------------------------------------------------------------------------
# Plot all the nodes in the given data base, and then color the control/ghost
# nodes of the given boundary condition independently.
#-------------------------------------------------------------------------------
def plotBoundaryNodes(dataBase, boundary):

    # First build one set of position pairs for all of the nodes in the
    # data base.
    positions = []
    for nodeList in dataBase.nodeLists:
        for r in list(nodeList.positions())[:nodeList.numInternalNodes]:
            positions.append((r.x, r.y))

    # Now build a list of the control node positions from the boundary
    # condition.
    controlPositions = []
    for nodeList in dataBase.nodeLists:
        controlNodes = boundary.controlNodes(nodeList)
        for nodeID in controlNodes:
            r = nodeList.positions()[nodeID]
            controlPositions.append((r.x, r.y))

    # Now build a list of the ghost node positions from the boundary
    # condition.
    ghostPositions = []
    for nodeList in dataBase.nodeLists:
        ghostNodes = boundary.ghostNodes(nodeList)
        for nodeID in ghostNodes:
            r = nodeList.positions()[nodeID]
            ghostPositions.append((r.x, r.y))

    # Finally we can plot these various sets of nodes.
    plot = plotXYTuples([positions, controlPositions, ghostPositions])
    return plot

#-------------------------------------------------------------------------------
# Plot the given sequences of (x,y) pairs, each with a distinct color.
#  [ [(x0,y0), (x1,y1), ...],
#    [(x0,y0), (x1,y1), ...],
#    .
#    .
#    .
#    [(x0,y0), (x1,y1), ...] ]
#-------------------------------------------------------------------------------
def plotXYTuples(listOfXYTuples):

    # Find the (min,max) of X and Y for all sets.
    xmin, ymin, xmax, ymax = findPairMinMax(listOfXYTuples[0])
    for seq in listOfXYTuples[1:]:
        xmin0, ymin0, xmax0, ymax0 = findPairMinMax(seq)
        xmin = min(xmin, xmin0)
        ymin = min(ymin, ymin0)
        xmax = max(xmax, xmax0)
        ymax = max(ymax, ymax0)

    # Create our plot result.
    plot = plt.figure()
    plot.axes().set_aspect("equal", "datalim")

    # Loop over the list of sequences of positions.
    color = iter(pltcm.rainbow(np.linspace(0,1,len(listOfXYTuples))))
    for seq in listOfXYTuples:
        c = next(color)

        # Build the local arrays of x and y.
        x = np.array([0.0]*len(seq))
        y = np.array([0.0]*len(seq))
        for i in range(len(seq)):
            x[i] = seq[i][0]
            y[i] = seq[i][1]

        # Plot this set of data.
        plot.plot(x, y, marker = "o", color = c)

    # That"s it, return the plot.
    return plot

#-------------------------------------------------------------------------------
# Find the (min, max) of a set of pairs.
#-------------------------------------------------------------------------------
def findPairMinMax(listOfPairs):
    minX, minY = 1e90, 1e90
    maxX, maxY = -1e90, -1e90
    for pair in listOfPairs:
        minX = min(minX, pair[0])
        minY = min(minY, pair[1])
        maxX = max(maxX, pair[0])
        maxY = max(maxY, pair[1])
    return minX, minY, maxX, maxY

#-------------------------------------------------------------------------------
# Plot the velocity field as a set of arrows.
# This is maintained here for backward compatibility, as a specialization of
# plotVectorField2d.
#-------------------------------------------------------------------------------
def plotVelocityField2d(dataBase,
                        plotGhosts = False,
                        velMultiplier = 1.0,
                        colorNodeLists = False,
                        colorDomains = False,
                        title = ""):

    return plotVectorField2d(dataBase,
                             dataBase.globalVelocity,
                             plotGhosts,
                             velMultiplier,
                             colorNodeLists,
                             colorDomains,
                             title)

#-------------------------------------------------------------------------------
# Plot the node spacing in 1D.
#-------------------------------------------------------------------------------
def plotNodeSpacing1d(dataBase):
    pos = dataBase.globalPosition
    xvals = []
    for ifield in range(len(pos)):
        xvals += [pos[ifield][i].x for i in range(pos[ifield].numInternalElements)]
    xvals = mpi.allreduce(xvals, mpi.SUM)
    xvals.sort()
    deltas = [xvals[i+1] - xvals[i] for i in range(len(xvals) - 1)] + [xvals[-1] - xvals[-2]]
    plot = newFigure()
    if plot:
        plot.plot(xvals, deltas)
    return plot

#-------------------------------------------------------------------------------
# Plot an arbitrary vector field as a set of arrows.
#-------------------------------------------------------------------------------
def plotVectorField2d(dataBase, fieldList,
                      plotGhosts = False,
                      vectorMultiplier = 1.0,
                      colorNodeLists = False,
                      colorDomains = False,
                      title = ""):

    assert colorNodeLists + colorDomains <= 1

    # Gather the node positions and vectors across all domains.
    # Loop over all the NodeLists.
    localNumNodes = []
    xNodes = []
    yNodes = []
    vxNodes = []
    vyNodes = []
    for i in range(dataBase.numNodeLists):
        nodeList = dataBase.nodeLists[i]
        assert i < fieldList.numFields
        vectorField = fieldList[i]
        if plotGhosts:
            n = nodeList.numNodes
        else:
            n = nodeList.numInternalNodes
        localNumNodes.append(n)
        xNodes += np.array([x.x for x in nodeList.positions()[:n]])
        yNodes += np.array([x.y for x in nodeList.positions()[:n]])
        vxNodes += np.array([x.x for x in vectorField[:n]])
        vyNodes += np.array([x.y for x in vectorField[:n]])
    assert len(xNodes) == len(yNodes) == len(vxNodes) == len(vyNodes)
    
    numDomainNodes = [len(xNodes)]
    numNodesPerDomain = mpi.gather(numDomainNodes)
    globalNumNodes = mpi.gather(localNumNodes)
    globalXNodes = mpi.gather(xNodes)
    globalYNodes = mpi.gather(yNodes)
    globalVxNodes = mpi.gather(vxNodes)
    globalVyNodes = mpi.gather(vyNodes)

    plot = newFigure()
    if mpi.rank == 0:
        plot.axes().set_aspect("equal", "datalim")
        plot.title(title)
        color = iter(pltcm.rainbow(np.linspace(0,1,len(xlist))))

        if colorDomains:
            cumulativeN = 0
            for domain in range(len(numNodesPerDomain)):
                c = next(color)
                n = numNodesPerDomain[domain]
                x = np.array(globalXNodes[cumulativeN:cumulativeN + n])
                y = np.array(globalYNodes[cumulativeN:cumulativeN + n])
                vx = np.array(globalVxNodes[cumulativeN:cumulativeN + n])
                vy = np.array(globalVyNodes[cumulativeN:cumulativeN + n])
                cumulativeN += n
                plot.quiver(x, y, vx, vy, color = c)

        elif colorNodeLists:
            cumulativeN = 0
            for i in range(len(globalNumNodes)):
                c = next(color)
                n = globalNumNodes[i]
                if n > 0:
                    iNodeList = i % dataBase.numNodeLists
                    x = np.array(globalXNodes[cumulativeN:cumulativeN + n])
                    y = np.array(globalYNodes[cumulativeN:cumulativeN + n])
                    vx = p.array(globalVxNodes[cumulativeN:cumulativeN + n])
                    vy = np.array(globalVyNodes[cumulativeN:cumulativeN + n])
                    cumulativeN += n
                    plot.quiver(x, y, vx, vy, color = c)

        else:
            x = np.array(globalXNodes)
            y = np.array(globalYNodes)
            vx = np.array(globalVxNodes)
            vy = np.array(globalVyNodes)
            plot.quiver(x, y, vx, vy)

    return plot

#-------------------------------------------------------------------------------
# Generate a regularly spaced sampling of the given FieldList
# The answer is returned in a 2-D numpy array.
#-------------------------------------------------------------------------------
def gridSample(fieldList,
               zFunction = "%s",
               nx = 100,
               ny = 100,
               xmin = None,
               xmax = None,
               ymin = None,
               ymax = None):

    assert nx > 0 and ny > 0

    # Set up our return value array.
    xValues = np.array([[0.0]*nx]*ny)
    yValues = np.array([[0.0]*nx]*ny)
    zValues = np.array([[0.0]*nx]*ny)

    # Gather the fieldList info across all processors to process 0.
    localNumNodes = []
    localX = []
    localY = []
    for ifield in range(fieldList.numFields):
        field = fieldList[ifield]
        n = field.nodeList().numNodes
        localNumNodes.append(n)
        for r in field.nodeList().positions():
            localX.append(r.x)
            localY.append(r.y)
    globalNumNodes = mpi.gather(localNumNodes)
    globalX = mpi.gather(localX)
    globalY = mpi.gather(localY)

    # If the user did not specify the sampling volume, then find the min and
    # max node positions.
    if xmin == None:
        xmin = min(localX)
    if ymin == None:
        ymin = min(localY)
    if xmax == None:
        xmax = max(localX)
    if ymax == None:
        ymax = max(localY)
    xmin = mpi.allreduce(xmin, mpi.MIN)
    ymin = mpi.allreduce(ymin, mpi.MIN)
    xmax = mpi.allreduce(xmax, mpi.MAX)
    ymax = mpi.allreduce(ymax, mpi.MAX)

    assert xmax > xmin
    assert ymax > ymin

    # Figure out the sizes of the bins we're going to be sampling in
    dx = (xmax - xmin)/nx
    dy = (ymax - ymin)/ny

    # Loop over all the grid sampling positions, and figure out this processors
    # contribution.
    for iy in range(ny):
        for ix in range(nx):
            xValues[iy][ix] = xmin + (ix + 0.5)*dx
            yValues[iy][ix] = ymin + (iy + 0.5)*dy
            r = Vector2d(xValues[iy][ix], yValues[iy][ix])
            z = fieldList.sample(r)
            localZ = eval(zFunction % "z")
            globalZ = mpi.reduce(localZ, mpi.SUM)
            if mpi.rank == 0:
                print("%i %i %i %s %g %g" % (mpi.rank, ix, iy, r, z, localZ))
                print("%i %g" % (mpi.rank, globalZ))
                zValues[iy][ix] = globalZ

    return xValues, yValues, zValues

#-------------------------------------------------------------------------------
# Plot the energy history of the given conservation object.
#-------------------------------------------------------------------------------
def plotEHistory(conserve):
    plot = newFigure()
    if mpi.rank == 0:
        t = conserve.timeHistory
        E = conserve.EHistory
        KE = conserve.KEHistory
        TE = conserve.TEHistory
        UE = conserve.EEHistory
        plot.plot(t, E, "k-", label="Total energy")
        plot.plot(t, KE, "b-", label="Kinetic energy")
        plot.plot(t, TE, "r-", label="Thermal energy")
        plot.plot(t, UE, "g-", label="Potential energy")
        plot.axes.legend()
    return plot

#-------------------------------------------------------------------------------
# Plot the linear momentum history of the given conservation object.
#-------------------------------------------------------------------------------
def plotpmomHistory(conserve):
    plot = newFigure()
    if mpi.rank == 0:
        t = conserve.timeHistory
        p = conserve.pmomHistory
        px = [x.x for x in p]
        py = [x.y for x in p]
        pz = [x.z for x in p]
        pmag = [x.magnitude() for x in p]
        plot.plot(t, px, "b-", label = "x momentum")
        plot.plot(t, py, "r-", label = "y momentum")
        plot.plot(t, pz, "g-", label = "z momentum")
        plot.plot(t, pmag, "b-", label = "total momentum")
    return plot

#-------------------------------------------------------------------------------
# Plot a surface
#-------------------------------------------------------------------------------
def plotSurface(x,   # 2D numpy array with x-coordinates        : shape (nx,ny)
                y,   # 2D numpy array with y-coordinates        : shape (nx,ny)
                z,   # 2D numpy array with z values for surface : shape (nx,ny)
                cmap = pltcm.coolwarm,   # Colormap
                xlabel = None,
                ylabel = None,
                zlabel = None,
                title = None):
    fig, ax = plt.subplots(subplot_kw = {"projection" : "3d"})
    surf = ax.plot_surface(x, y, z, cmap = cmap)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_zlabel(zlabel)
    fig.colorbar(surf)
    plt.title(title)
    return fig, ax, surf

#-------------------------------------------------------------------------------
# Plot a QuadraticInterpolator
#-------------------------------------------------------------------------------
def plotInterpolator(interp,
                     n = None,
                     plot = None,
                     plotstyle = "r-",
                     label = None,
                     xlabel = None,
                     ylabel = None,
                     title = None):
    x0, x1 = interp.xmin, interp.xmax
    if n is None:
        n = 2 * interp.size
    if plot is None:
        plot = newFigure()
    xvals = np.linspace(x0, x1, n)
    yvals = np.array([interp(x) for x in xvals])
    plot.plot(xvals, yvals, plotstyle, label=label)
    plot.set_xlabel(xlabel)
    plot.set_ylabel(ylabel)
    plot.set_title(title)
    return plot

#-------------------------------------------------------------------------------
# Plot a table kernel
#-------------------------------------------------------------------------------
def plotTableKernel(WT, nPerh):
    plots = [plotInterpolator(interp = x,
                              xlabel = xlab,
                              ylabel = ylab,
                              title = ylab) for x, xlab, ylab in [(WT.Winterpolator,      r"$\eta$",   r"$W(\eta)$"),
                                                                  (WT.gradWinterpolator,  r"$\eta$",   r"$\partial_\eta W(\eta)$"),
                                                                  (WT.grad2Winterpolator, r"$\eta$",   r"$\partial^2_\eta W(\eta)$"),
                                                                  (WT.nPerhInterpolator,  r"$\sum W$",  r"n per h($\sum W$)"),
                                                                  (WT.WsumInterpolator,   r"n per h",  r"$\sum W$")]]

    x0, x1 = 0.0, WT.kernelExtent
    xvals = np.linspace(x0, x1, 100)
    yvals = np.array([WT.kernelValueSPH(x) for x in xvals])
    plotSPH = newFigure()
    plotSPH.plot(xvals, yvals, "r-", label=None)
    plotSPH.set_xlabel(r"$\eta$")
    plotSPH.set_ylabel(r"$W_{SPH}(\eta)$")
    plotSPH.set_title(r"$W(\eta)$ for SPH h lookup")

    yvals = np.array([WT.kernelValueASPH(x, nPerh) for x in xvals])
    plotASPH = newFigure()
    plotASPH.plot(xvals, yvals, "r-", label=None)
    plotASPH.set_xlabel(r"$\eta$")
    plotASPH.set_ylabel(r"$W_{ASPH}(\eta)$")
    plotASPH.set_title(r"$W(\eta)$ for ASPH h lookup with $n_h="+str(nPerh)+"$")

    plots += [plotSPH, plotASPH]

    return plots

#-------------------------------------------------------------------------------
# Plot a polygon.
#-------------------------------------------------------------------------------
def plotPolygon(polygon,
                plotVertices = True,
                plotFacets = True,
                plotNormals = False,
                plotCentroid = False,
                plot = None,
                persist = False,
                plotLabels = True):
    mppoly = patches.Polygon(np.array([[v.x, v.y] for v in polygon.vertices]), False)

    if plot is None:
        plot = newFigure()
    plot.add_patch(mppoly)
    return

    # px = []
    # py = []
    # for v in polygon.vertices():
    #     px.append(v.x)
    #     py.append(v.y)
    # fx = []
    # fy = []
    # fdx = []
    # fdy = []
    # nx = []
    # ny = []
    # ndx = []
    # ndy = []
    # for f in polygon.facets():
    #     dr = f.point2 - f.point1
    #     hdr = dr/2.0
    #     fx.append(f.point1.x)
    #     fy.append(f.point1.y)
    #     fdx.append(dr.x)
    #     fdy.append(dr.y)
    #     nx.append(fx[-1] + hdr.x)
    #     ny.append(fy[-1] + hdr.y)
    #     ndx.append(f.normal.x)
    #     ndy.append(f.normal.y)
    # if plot is None:
    #     plot = generateNewGnuPlot(persist)
    # if plotLabels:
    #     vlabel, flabel, nlabel = "Vertices", "Facets", "Normals"
    # else:
    #     vlabel, flabel, nlabel = None, None, None
    # dataPoints = Gnuplot.Data(px, py,
    #                           with_ = "points pt 1 ps 2",
    #                           title = vlabel,
    #                           inline = True)
    # dataFacets = Gnuplot.Data(fx, fy, fdx, fdy,
    #                           with_ = "vectors",
    #                           title = flabel,
    #                           inline = True)
    # dataNormals = Gnuplot.Data(nx, ny, ndx, ndy,
    #                            with_ = "vectors",
    #                            title = nlabel,
    #                            inline = True)
    # if plotVertices:
    #     plot.replot(dataPoints)

    # if plotFacets:
    #     plot.replot(dataFacets)

    # if plotNormals:
    #     plot.replot(dataNormals)

    # if plotCentroid:
    #     c = polygon.centroid()
    #     dataCentroid = Gnuplot.Data([c.x], [c.y],
    #                                 with_ = "points pt 2 ps 2",
    #                                 title = "Centroid",
    #                                 inline = True)
    #     plot.replot(dataCentroid)

    # SpheralGnuPlotCache.extend([dataPoints, dataFacets, dataNormals, plot])

    # return plot

# #-------------------------------------------------------------------------------
# # Plot a PolygonalMesh
# #-------------------------------------------------------------------------------
# def plotPolygonalMesh(mesh,
#                       persist = False):
#     polylocal = []
#     for izone in xrange(mesh.numZones):
#         zone = mesh.zone(izone)
#         polylocal.append([mesh.node(i).position() for i in zone.nodeIDs])
#         polylocal[-1].append(polylocal[-1][0])
#     assert len(polylocal) == mesh.numZones

#     p = generateNewGnuPlot(persist)
#     for sendProc in xrange(mpi.procs):
#         polys = mpi.bcast(polylocal, root=sendProc)
#         for poly in polys:
#             p.replot(Gnuplot.Data([x.x for x in poly], [x.y for x in poly],
#                                   with_ = "lines lt %i lw 2" % 1,
#                                   title = None,
#                                   inline = True))
#     return p
