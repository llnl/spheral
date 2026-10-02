#-------------------------------------------------------------------------------
# VoronoiCells
#-------------------------------------------------------------------------------
from PYB11Generator import *
from VolumeUpdate import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11dynamic_attr
class VoronoiCells(VolumeUpdate):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using Vector = typename %(Dimension)s::Vector;
    using Tensor = typename %(Dimension)s::Tensor;
    using SymTensor = typename %(Dimension)s::SymTensor;
    using FacetedVolume = typename %(Dimension)s::FacetedVolume;
    using TimeStepType = typename Physics<%(Dimension)s>::TimeStepType;
    using ResidualType = typename Physics<%(Dimension)s>::ResidualType;
    using VolumeRequirements = typename Physics<%(Dimension)s>::VolumeRequirements;
    using RKRequirements = typename Physics<%(Dimension)s>::RKRequirements;
    using ConnectivityRequirements = typename Physics<%(Dimension)s>::ConnectivityRequirements;
"""
    
    def pyinit(self,
               volumeType = ("const VolumeType", "VolumeType::VoronoiVolume"),
               W = "const TableKernel<%(Dimension)s>&",
               facetedBoundaries = ("const std::vector<FacetedVolume>&", "std::vector<FacetedVolume>()"),
               facetedHoles = ("const std::vector<std::vector<FacetedVolume>>&", "std::vector<std::vector<FacetedVolume>>()"),
               updateInStep = ("const bool", "true"),
               updateInFinalize = ("const bool", "false")):
        "VoronoiCells constructor"
        return

    #...........................................................................
    # Methods unique to VoronoiCells
    @PYB11virtual
    def addFacetedBoundary(bound = "const FacetedVolume&",
                           holes = ("const std::vector<FacetedVolume>&", "std::vector<FacetedVolume>()")):
        "Add a faceted boundary (optionally with holes)"
        return "void"

    #...........................................................................
    # Properties (Voronoi-specific; volume/updateInStep/etc inherited from VolumeUpdate)
    kernelExtent =      PYB11property("Scalar",                                                     "kernelExtent", doc="The kernel extent in eta")
    weight =            PYB11property("const FieldList<%(Dimension)s, Scalar>&",                    "weight",            returnpolicy="reference_internal")
    surfacePoint =      PYB11property("const FieldList<%(Dimension)s, int>&",                       "surfacePoint",      returnpolicy="reference_internal")
    etaVoidPoints =     PYB11property("const FieldList<%(Dimension)s, std::vector<Vector>>&",       "etaVoidPoints",     returnpolicy="reference_internal")
    cells =             PYB11property("const FieldList<%(Dimension)s, FacetedVolume>&",             "cells",             returnpolicy="reference_internal")
    cellFaceFlags =     PYB11property("const FieldList<%(Dimension)s, std::vector<CellFaceFlag>>&", "cellFaceFlags",     returnpolicy="reference_internal")
    deltaCentroid =     PYB11property("const FieldList<%(Dimension)s, Vector>&",                    "deltaCentroid",     returnpolicy="reference_internal")
    facetedBoundaries = PYB11property("const std::vector<FacetedVolume>&",                          "facetedBoundaries", returnpolicy="reference_internal")
    facetedHoles =      PYB11property("const std::vector<std::vector<FacetedVolume>>&",             "facetedHoles",      returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, VoronoiCells)
