"""
Spheral NodeGenerators module.

Algorithms and methods for creating initial conditions in Spheral.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from WeightingFunctor import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"NodeGenerators/generateCylDistributionFromRZ.hh"',
                  '"NodeGenerators/fillFacetedVolume.hh"',
                  '"NodeGenerators/relaxNodeDistribution.hh"',
                  '"NodeGenerators/readSiloPolyMesh.hh"',
                  '"NodeGenerators/centroidalRelaxNodesImpl.hh"',
                  '"NodeGenerators/compactFacetedVolumes.hh"',
                  '"NodeGenerators/chooseRandomNonoverlappingCenter.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate types and add dimension dependent functions.
#-------------------------------------------------------------------------------
if 3 in dims:
    def fillFacetedVolume(outerBoundary = "const Dim<3>::FacetedVolume&",
                          n1d = "const unsigned",
                          domain = "const unsigned",
                          numDomains = "const unsigned"):
        "Fill an outer bounding volume (specify x number of points)."
        return "std::vector<Dim<3>::Vector>"

    def fillFacetedVolume2(outerBoundary = "const Dim<3>::FacetedVolume&",
                           dx = "const double",
                           domain = "const unsigned",
                           numDomains = "const unsigned"):
        "Fill an outer bounding volume (dx specified)."
        return "std::vector<Dim<3>::Vector>"

    def fillFacetedVolume3(innerBoundary = "const Dim<3>::FacetedVolume&",
                           outerBoundary = "const Dim<3>::FacetedVolume&",
                           n1d = "const unsigned",
                           domain = "const unsigned",
                           numDomains = "const unsigned"):
        "Fill between an inner and outer boundary (specify x number of points)."
        return "std::vector<Dim<3>::Vector>"

    def fillFacetedVolume10(outerBoundary = "const Dim<3>::FacetedVolume&",
                            innerBoundary = "const Dim<3>::FacetedVolume&",
                            dx = "const double",
                            domain = "const unsigned",
                            numDomains = "const unsigned"):
        "Fill between an inner and outer boundary (dx specified)."
        return "std::vector<Dim<3>::Vector>"

    def readSiloPolyMesh(fileName = "const std::string&",
                         meshName = "const std::string&",
                         positions = "std::vector<Dim<3>::Vector>&",
                         volumes = "std::vector<double>&",
                         H = "std::vector<Dim<3>::SymTensor>&"):
        """Helper method for the SiloPolyMeshGenerator which reads a polyhedral mesh 
from a silo file and returns the geometry."""
        return "void"

#...............................................................................
if 2 in dims and 3 in dims:
    def generateCylDistributionFromRZ(x = "std::vector<double>&",
                                      y = "std::vector<double>&",
                                      z = "std::vector<double>&",
                                      m = "std::vector<double>&",
                                      H = "std::vector<Dim<3>::SymTensor>&",
                                      globalIDs = "std::vector<int>&",
                                      extraFields = "std::vector<std::vector<double> >&",
                                      nNodePerh = "const double",
                                      kernelExtent = "const double",
                                      phi = "const double",
                                      procID = "const int",
                                      nProcs = "const int"):
        """Helper method for the GenerateCylindricalNodeDistribution3d node generator
to generate the spun node distribution."""
        return "void"

#...............................................................................
@PYB11template("Dimension")
def centroidalRelaxNodesImpl(db = "DataBase<%(Dimension)s>&",
                             volumeBoundaries = "const std::vector<typename %(Dimension)s::FacetedVolume>&",
                             holes = "const std::vector<std::vector<typename %(Dimension)s::FacetedVolume> >&",
                             W = "const TableKernel<%(Dimension)s>&",
                             rhofunc = "const PythonBoundFunctors::SpheralFunctor<typename %(Dimension)s::Vector, double>&",
                             gradrhofunc = "const PythonBoundFunctors::SpheralFunctor<typename %(Dimension)s::Vector, typename %(Dimension)s::Vector>&",
                             rhoConst = "const bool",
                             useGradRhoFunc = "const bool",
                             boundaries = "std::vector<Boundary<%(Dimension)s>*>&",
                             maxIterations = "const unsigned",
                             maxFracTol = "const double",
                             avgFracTol = "const double",
                             correctionOrder = "const RKOrder",
                             centroidFrac = "const double",
                             vol = "FieldList<%(Dimension)s, double>&",
                             surfacePoint = "FieldList<%(Dimension)s, int>&",
                             cells = "FieldList<%(Dimension)s, typename %(Dimension)s::FacetedVolume>&"):
    "Implement Lloyd's algorithm for centroidal relaxation of fluid points."
    return "unsigned"

for ndim in dims:
    exec('''
centroidalRelaxNodesImpl%(ndim)id = PYB11TemplateFunction(centroidalRelaxNodesImpl, template_parameters="%(Dimension)s", pyname="centroidalRelaxNodesImpl")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

#...............................................................................
@PYB11template("Dimension")
def relaxNodeDistribution(dataBase = "DataBase<%(Dimension)s>&",
                          boundary = "const typename %(Dimension)s::FacetedVolume&",
                          boundaries = "const std::vector<Boundary<%(Dimension)s>*>&",
                          W = "const TableKernel<%(Dimension)s>&",
                          weightingFunctor = "const WeightingFunctor<%(Dimension)s>&",
                          massDensityFunctor = "const WeightingFunctor<%(Dimension)s>&",
                          targetMass = "const double",
                          maxIterations = "const int",
                          tolerance = "const double"):
    """Centroidally relax a node distribution in a boundary.
Optionally the user can specify a weighting function for the nodes."""
    return "void"

@PYB11template("Dimension")
def compactFacetedVolumes(shapes = "std::vector<typename %(Dimension)s::FacetedVolume>&",
                          centers = "std::vector<typename %(Dimension)s::Vector>&",
                          flags = "std::vector<int>&",
                          surface = "const typename %(Dimension)s::FacetedVolume&",
                          depthmax = "const double",
                          surfaceIterations = "const unsigned",
                          maxIterations = "const unsigned",
                          dispfrac = "const double",
                          maxoverlapfrac = "const double"):
    "Push FacetedVolume shapes together inside a surface, but excluding mutual overlap."
    return "unsigned"

@PYB11template("Dimension")
def chooseRandomNonoverlappingCenter(result = "typename %(Dimension)s::Vector&",
                                     trialShape = "const typename %(Dimension)s::FacetedVolume&",
                                     boundary = "const typename %(Dimension)s::FacetedVolume&",
                                     existingShapes = "const std::vector<typename %(Dimension)s::FacetedVolume>&",
                                     maxTries = "const unsigned"):
    return "unsigned"

subdims = [x for x in (2, 3) if x in dims]
for ndim in subdims:
    exec('''
WeightingFunctor%(ndim)id = PYB11TemplateClass(WeightingFunctor, template_parameters="%(Dimension)s")
relaxNodeDistribution%(ndim)id = PYB11TemplateFunction(relaxNodeDistribution, template_parameters="%(Dimension)s", pyname="relaxNodeDistribution")
compactFacetedVolumes%(ndim)id = PYB11TemplateFunction(compactFacetedVolumes, template_parameters="%(Dimension)s", pyname="compactFacetedVolumes")
chooseRandomNonoverlappingCenter%(ndim)id = PYB11TemplateFunction(chooseRandomNonoverlappingCenter, template_parameters="%(Dimension)s", pyname="chooseRandomNonoverlappingCenter")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
