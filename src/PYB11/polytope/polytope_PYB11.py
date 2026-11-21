"""
polytope module.

Rewrap the polytope interface for Python with pybind11.  This should
be superseeded by natively redoing polytope's Python interface with 
PYB11 at some point.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

PYB11includes = ['"polytope/Tessellation.hh"',
                 '"Mesh/copy2polytope.hh"',
                 '"Geometry/GeomPlane.hh"',
                 '<sstream>']

PYB11namespaces = ["Spheral",
                   "polytope"]

PYB11opaque = ["std::vector<char>",
               "std::vector<unsigned>",
               "std::vector<uint64_t>",
               "std::vector<int>",
               "std::vector<float>",
               "std::vector<double>",
               "std::vector<std::string>",

               "std::vector<std::vector<char>>",
               "std::vector<std::vector<unsigned>>",
               "std::vector<std::vector<uint64_t>>",
               "std::vector<std::vector<int>>",
               "std::vector<std::vector<float>>",
               "std::vector<std::vector<double>>",
               "std::vector<std::vector<std::string>>"]

@PYB11template("Dimension", "nDim")
def copy2polytope(cells = "const FieldList<%(Dimension)s, %(Dimension)s::FacetedVolume>&",
                  mesh = "polytope::Tessellation<%(nDim)s, double>&"):
    "Copy a FieldList of polygons/polyhedra to a polytope Tessellation."
    return "void"

# Some methods are only valid in 2D and 3D.
for ndim in (x for x in dims if x in (2, 3)):
    exec('''
copy2polytope%(ndim)id = PYB11TemplateFunction(copy2polytope, template_parameters=("Dim<%(ndim)i>", "%(ndim)i"), pyname="copy2polytope")
''' % {"ndim"   : ndim})
