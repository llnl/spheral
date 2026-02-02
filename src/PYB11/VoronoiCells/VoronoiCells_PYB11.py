"""
Spheral VoronoiCells module.

Provides VoronoiCells Spheral package
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from VoronoiCells import *
from SubPointPressureHourglassControl import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"VoronoiCells/VoronoiCells.hh"',
                  '"VoronoiCells/SubPointPressureHourglassControl.hh"',
                  '"VoronoiCells/IncrementVoronoiCells.hh"',
                  '"VoronoiCells/computeVoronoiVolume.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"FileIO/FileIO.hh"',
                  '"Boundary/Boundary.hh"',
                  '<sstream>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeVoronoiVolume(position = "const FieldList<%(Dimension)s, %(Dimension)s::Vector>&",
                         H = "const FieldList<%(Dimension)s, %(Dimension)s::SymTensor>&",
                         connectivityMap = "const ConnectivityMap<%(Dimension)s >&",
                         damage = "const FieldList<%(Dimension)s, %(Dimension)s::SymTensor>&",
                         facetedBoundaries = "const std::vector<%(Dimension)s::FacetedVolume>&",
                         holes = "const std::vector<std::vector<%(Dimension)s::FacetedVolume> >&",
                         boundaries = "const std::vector<Boundary<%(Dimension)s>*>&",
                         weight = "const FieldList<%(Dimension)s, %(Dimension)s::Scalar>&",
                         surfacePoint = "FieldList<%(Dimension)s, int>&",
                         vol = "FieldList<%(Dimension)s, %(Dimension)s::Scalar>&",
                         deltaMedian = "FieldList<%(Dimension)s, %(Dimension)s::Vector>&",
                         etaVoidPoints = "FieldList<%(Dimension)s, std::vector<%(Dimension)s::Vector>>&",
                         cells = "FieldList<%(Dimension)s, %(Dimension)s::FacetedVolume>&",
                         cellFaceFlags = "FieldList<%(Dimension)s, std::vector<CellFaceFlag>>&"):
    "Compute the volume per point based on the Voronoi tessellation-like algorithm."
    return "void"

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    Dimension = "Dim<{}>".format(ndim)
    exec(f'''
computeVoronoiVolume{ndim}d = PYB11TemplateFunction(computeVoronoiVolume, template_parameters="{Dimension}", pyname="computeVoronoiVolume")

VoronoiCells{ndim}d = PYB11TemplateClass(VoronoiCells, template_parameters="{Dimension}")
SubPointPressureHourglassControl{ndim}d = PYB11TemplateClass(SubPointPressureHourglassControl, template_parameters="{Dimension}")
''')

    # % {ndim      : ndim,
    #    Dimension : "Dim<{}>".format(ndim)})
