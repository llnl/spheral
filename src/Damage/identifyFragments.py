import Spheral
import mpi

import sys
import time
from SpheralTestUtilities import fuzzyEqual

#-------------------------------------------------------------------------------
# Helper method to get the unique items in a list.
#-------------------------------------------------------------------------------
def uniqueElements(container):
    resultDict = {}
    for x in container:
        resultDict[x] = 0
    return list(resultDict.keys())

#-------------------------------------------------------------------------------
# Identify the set of distinct fragments in the given SolidNodeList.
#
# Returns a field, filled with integer indicies for each node identifying the
# fragment it is associated with.
#-------------------------------------------------------------------------------
def identifyFragments(nodeList,
                      linkRadius,
                      damageField,
                      damageThreshold,
                      assignDustToFragments,
                      density = None,
                      densityThreshold = None):

    # Preconditions.
    assert linkRadius > 0.0

    startTime = time.time()

    # Decide what our dimensionality is.
    if isinstance(nodeList, Spheral.NodeList1d):
        FieldConstructor = Spheral.IntField1d
        computeFragMethod = Spheral.computeFragmentField1d
    elif isinstance(nodeList, Spheral.NodeList2d):
        FieldConstructor = Spheral.IntField2d
        computeFragMethod = Spheral.computeFragmentField2d
    elif isinstance(nodeList, Spheral.NodeList3d):
        FieldConstructor = Spheral.IntField3d
        computeFragMethod = Spheral.computeFragmentField3d
    else:
        raise "identifyFragments ERROR: What the heck is %s!  I expected a NodeList." % str(nodeList)

    # Check for any mass density criteria
    if density is None:
        density = nodeList.massDensity()
    if densityThreshold is None:
        densityThreshold = 0.0

    # Eliminate any ghost nodes.
    nodeList.numGhostNodes = 0

    # The C++ method does the heavy lifting.
    # Compute the fragments on the undamaged NodeList.
    mask = nodeList.mask();
    result = computeFragMethod(nodeList, linkRadius, density, damageField, mask, densityThreshold, damageThreshold, assignDustToFragments)
    numFragments = mpi.allreduce(max(list(result.internalValues()) + [0]), mpi.MAX)
    numFragments += 1

    # To avoid confusion with the stored Field in SolidNodeList, we rename the result
    # of this method.
    result.name = "FRAGMENT INDEX"

    stopTime = time.time()
    print("identifyFragments: identified %i fragments." % numFragments)
    print("                   Required %g seconds." % (stopTime - startTime))

    return result

#-------------------------------------------------------------------------------
# Determine some interesting properties for a collection of fragments.
#
# Result is a list (indexed by fragment ID) of dictionaries for each fragment,
# consisting of:
#   mass
#   num nodes
#   position (center of mass)
#   velocity (center of mass)
#   volume
#   average mass density
#   average thermal energy
#   average pressure
#   shape ellipsoid tensor (all nodes)
#   shape ellipsoid tensor (damage thresholded)
#-------------------------------------------------------------------------------
def fragmentProperties(nodeList,
                       fragField,
                       strain = None):

    startTime = time.time()

    # Decide what our dimensionality is.
    if isinstance(nodeList, Spheral.NodeList1d):
        Vector = Spheral.Vector1d
        SymTensor = Spheral.SymTensor1d
        ScalarField = Spheral.ScalarField1d
    elif isinstance(nodeList, Spheral.NodeList2d):
        Vector = Spheral.Vector2d
        SymTensor = Spheral.SymTensor2d
        ScalarField = Spheral.ScalarField2d
    elif isinstance(nodeList, Spheral.NodeList3d):
        Vector = Spheral.Vector3d
        SymTensor = Spheral.SymTensor3d
        ScalarField = Spheral.ScalarField3d
    else:
        raise "fragmentProperties ERROR: What the heck is %s!  I expected a NodeList." % str(nodeList)

    # Are we compiling stats on the strain?
    usingStrain = not (strain is None)

    # Determine how many fragments there are.
    numFragments = mpi.allreduce(max(list(fragField.internalValues()) + [0]), mpi.MAX) + 1
    assert numFragments >= 1

    # Prepare the result.
    result = {}
    for i in range(numFragments):
        result[i] = {"mass": 0.0,
                     "num nodes": 0,
                     "position": Vector(),
                     "velocity": Vector(),
                     "volume": 0.0,
                     "mass density": 0.0,
                     "thermal energy": 0.0,
                     "pressure": 0.0,
                     "shape tensor": SymTensor(),
                     "shape eigen": Vector(),
                     }
        if usingStrain:
            result[i]["strain (min)"] = 0.0
            result[i]["strain (max)"] = 0.0
            result[i]["strain (vol)"] = 0.0
    assert len(result) == numFragments

    # Grab the state fields.
    mass = nodeList.mass()
    position = nodeList.positions()
    velocity = nodeList.velocity()
    rho = nodeList.massDensity()
    u = nodeList.specificThermalEnergy()
    P = ScalarField("pressure", nodeList)
    nodeList.pressure(P)

    # Now iterate over the nodes and accumulate the local result.
    for i in range(nodeList.numInternalNodes):
        fragID = fragField[i]
        assert fragID < numFragments
        mi = mass[i]
        result[fragID]["mass"] += mi
        result[fragID]["num nodes"] += 1
        result[fragID]["position"] += position[i]*mi
        result[fragID]["velocity"] += velocity[i]*mi
        result[fragID]["volume"] += mi/rho[i]
        result[fragID]["mass density"] += mi*rho[i]
        result[fragID]["thermal energy"] += mi*u[i]
        result[fragID]["pressure"] += mi*P[i]
        if usingStrain:
            sev = strain[i].eigenValues()
            result[fragID]["strain (min)"] += mi*sev.minElement()
            result[fragID]["strain (max)"] += mi*sev.maxElement()
            result[fragID]["strain (vol)"] += mi*sev.sumElements()
    assert fuzzyEqual(sum([result[i]["mass"] for i in result]), sum(nodeList.mass().internalValues()))

    # Reduce for the global fragment properties.
    assert mpi.allreduce(len(result), mpi.SUM) == numFragments*mpi.procs
    for fragID in range(numFragments):
        mfrag = mpi.allreduce(result[fragID]["mass"], mpi.SUM)
        assert mfrag > 0.0
        result[fragID]["mass"] = mfrag
        result[fragID]["num nodes"] = mpi.allreduce(result[fragID]["num nodes"], mpi.SUM)
        result[fragID]["position"] = mpi.allreduce(result[fragID]["position"], mpi.SUM)/mfrag
        result[fragID]["velocity"] = mpi.allreduce(result[fragID]["velocity"], mpi.SUM)/mfrag
        result[fragID]["volume"] = mpi.allreduce(result[fragID]["volume"], mpi.SUM)
        result[fragID]["mass density"] = mpi.allreduce(result[fragID]["mass density"], mpi.SUM)/mfrag
        result[fragID]["thermal energy"] = mpi.allreduce(result[fragID]["thermal energy"], mpi.SUM)/mfrag
        result[fragID]["pressure"] = mpi.allreduce(result[fragID]["pressure"], mpi.SUM)/mfrag
        if usingStrain:
            result[fragID]["strain (min)"] = mpi.allreduce(result[fragID]["strain (min)"], mpi.SUM)/mfrag
            result[fragID]["strain (max)"] = mpi.allreduce(result[fragID]["strain (max)"], mpi.SUM)/mfrag
            result[fragID]["strain (vol)"] = mpi.allreduce(result[fragID]["strain (vol)"], mpi.SUM)/mfrag

    # Now that we have the center of mass for each fragment, we can evaluate the shape tensors.
    for i in range(nodeList.numInternalNodes):
        fragID = fragField[i]
        assert fragID < numFragments
        mi = mass[i]
        dr = position[i] - result[fragID]["position"]
        result[fragID]["shape tensor"] += dr.selfdyad() * mi

    # Reduce the global shapes.
    for fragID in range(numFragments):
        mfrag = result[fragID]["mass"]
        assert mfrag > 0.0
        result[fragID]["shape tensor"] = mpi.allreduce(result[fragID]["shape tensor"], mpi.SUM)/mfrag
        result[fragID]["shape eigen"] = result[fragID]["shape tensor"].eigenValues()

    # That's it.
    assert fuzzyEqual(sum([result[i]["mass"] for i in result]), mpi.allreduce(sum(nodeList.mass().internalValues()), mpi.SUM))

    stopTime = time.time()
    print("fragmentProperties:  Required %g seconds." % (stopTime - startTime))

    return result
