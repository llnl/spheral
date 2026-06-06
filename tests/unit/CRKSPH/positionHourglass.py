from Spheral1d import *

def positionHourglass(db, WT, boundaries, 
                      minfrac = 0.01,
                      maxfrac = 0.1):

    position = db.fluidPosition
    
    nodeLists = db.fluidNodeLists
    for nodes in nodeLists:
        nodes.numGhostNodes = 0
        nodes.neighbor().updateNodes()
    for bc in boundaries:
        bc.setAllGhostNodes(db)
        bc.finalizeGhostBoundary()
        for nodes in nodeLists:
            nodes.neighbor().updateNodes()

    db.updateConnectivityMap()
    cm = db.connectivityMap()
    position = db.fluidPosition
    mass = db.fluidMass
    rho = db.fluidMassDensity
    H = db.fluidHfield

    W0 = WT.kernelValue(0.0, 1.0)

    for bc in boundaries:
        bc.applyFieldListGhostBoundary(mass)
        bc.applyFieldListGhostBoundary(rho)
    for bc in boundaries:
        bc.finalizeGhostBoundary()

    result = 0
    for nodeListi in range(db.numFluidNodeLists):
        for i in range(nodeLists[nodeListi].numInternalNodes):
            connectivity = cm.connectivityForNode(nodeListi, i)
            for nodeListj in range(db.numFluidNodeLists):
                for j in connectivity[nodeListj]:
                    xji = position[nodeListj][j] - position[nodeListi][i]
                    xjihat = xji.unitVector()
                    delta = max(0.0, 0.5*(mass[nodeListi][i]/rho[nodeListi][i] + mass[nodeListj][j]/rho[nodeListj][j]) - xji.magnitude())
                    hi = 1.0/(H[nodeListi][i]*xjihat).magnitude()
                    mul = 1.0 # min(maxfrac, delta/hi)
                    etai = xji.magnitude()/hi
                    weight = WT.kernelValue(etai, 1.0)/W0
                    position[0][i] -= mul*weight*delta*xjihat

    for bc in boundaries:
        bc.setAllViolationNodes(db)
    for nodes in nodeLists:
        nodes.neighbor().updateNodes()

    return result

