from SpheralCompiledPackages import *

from spheralDimensions import spheralDimensions
dims = spheralDimensions()

# clean up the smooth w/ diffent normal calculations
#-------------------------------------------------------------------------------
# Convienent constructor for slide surfaces.
#-------------------------------------------------------------------------------
# RE LINE 20-21  
# if we have different smoothing lengths at a material interface w/ slip active, 
# the E0 error of the SPH kernel will cause the surface normals to have really 
# bad values two-or-so rows in from the interface. To get around this issue, 
# whenever we base our surface normal on its same-material neighbors we need
# an extra smoothing step to reorient normals a few rows back. So effectively
# we have an optionally default: 
# surfaceNormalMethod = DifferentMaterial
#   no smoothing
# surfaceNormalMethod = SameMaterial | AllMaterial | MassWeighted
#   smoothing
# this should be cleaned up once we have a better idea of what the best 
# approach is.
#-------------------------------------------------------------------------------
SlideSurfaceFactoryString = """
def makeSlideSurfaces%(dim)s(dataBase,
                             slideSurfaces=None):

        contactTypes = [0]*(dataBase.numNodeLists**2)

        if slideSurfaces:
            
            # create the map nodelist --> index
            nodeLists = dataBase.nodeLists
            nodeListMap = {}
            for i in range(dataBase.numNodeLists):        
                nodeListMap[nodeLists[i]]=i

            # table (2d->1d) this should be fixed later
            for slide in slideSurfaces:
                nodeListi = nodeListMap[slide[0]]
                nodeListj = nodeListMap[slide[1]]
                contactTypes[dataBase.numNodeLists*nodeListi+nodeListj]=1
                contactTypes[nodeListi+dataBase.numNodeLists*nodeListj]=1

        result = SlideSurface%(dim)s(dataBase,
                                     vector_of_int(contactTypes))

        return result
"""

for dim in dims:
    exec(SlideSurfaceFactoryString % {"dim": "%id" % dim}) 
