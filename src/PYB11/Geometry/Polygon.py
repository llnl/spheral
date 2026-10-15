#-------------------------------------------------------------------------------
# Polygon
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11cppname("GeomPolygon")
class Polygon:

    PYB11typedefs = """
    using Polygon = GeomPolygon;
    using Vector = GeomPolygon::Vector;
    using Tensor = GeomPolygon::Tensor;
    using Facet = GeomPolygon::Facet;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    def pyinit0(self,
                rhs = "const Polygon&"):
        "Copy constructor"

    def pyinit1(self,
                points = "const std::vector<Vector>&"):
        """Note this constructor constructs the convex hull of the given points,
meaning that the full set of points passed in may not appear in the vertices."""

    def pyinit2(self,
                points = "const std::vector<Vector>&",
                facetIndices = "const std::vector<std::vector<unsigned> >&"):
        "Construct with explicit vertices and facets"

    @PYB11implementation("[](py::list& points) { std::vector<Vector> vpoints; for (auto p: points) vpoints.push_back(p.cast<Vector>()); return new GeomPolygon(vpoints); }")
    def pyinit3(self,
                points = "py::list"):
        "Construct as the convex hull of a python list of points"

    #...........................................................................
    # Methods
    @PYB11const
    def contains(self,
                 point = "const Vector&",
                 countBoundary = ("const bool", "true"),
                 tol = ("const double", "1.0e-8")):
        "Test if the given point is internal to the polygon."
        return "bool"

    @PYB11const
    def convexContains(self,
                       point = "const Vector&",
                       countBoundary = ("const bool", "true"),
                       tol = ("const double", "1.0e-8")):
        "Test if the given point is internal to the polygon (assumes convexity)."
        return "bool"

    @PYB11const
    def intersect(self,
                  rhs = "const Polygon&"):
        "Test if we intersect another polygon."
        return "bool"

    @PYB11const
    def convexIntersect(self,
                        rhs = "const Polygon&"):
        "Test if we intersect another polygon (assumes convexity)"
        return "bool"

    @PYB11const
    @PYB11pycppname("intersect")
    def intersect1(self,
                   rhs = "const std::pair<Vector, Vector>&"):
        "Test if we intersect a box represented by a min/max pair of coordinates."
        return "bool"

    @PYB11const
    @PYB11pycppname("intersect")
    def intersect2(self,
                   s0 = "const Vector&",
                   s1 = "const Vector&",
                   countInterior = ("const bool", "true")):
        "Test if we intersect a line segment."
        return "bool"

    @PYB11const
    @PYB11implementation("[](const Polygon& self, const Vector& s0, const Vector& s1) { std::vector<unsigned> facetIDs; std::vector<Vector> intersections; self.intersections(s0, s1, facetIDs, intersections); return py::make_tuple(facetIDs, intersections); }")
    def intersections(self,
                      s0 = "const Vector&",
                      s1 = "const Vector&"):
        "Return the intersections of this polygon with a line segment denoted by it's end points."
        return "py::tuple"

    @PYB11const
    def edges(self):
        "Get the edges as integer (node) pairs."
        return "std::vector<std::pair<unsigned, unsigned> >"

    def reconstruct(self,
                    vertices = "const std::vector<Vector>&",
                    facetVertices = "const std::vector<std::vector<unsigned> >&"):
        """Reconstruct the internal data given a set of verticies and the vertex
indices that define the facets."""
        return "void"

    @PYB11const
    def closestFacet(self, p = "const Vector&"):
        "Find the facet closest to the given point."
        return "unsigned"

    @PYB11const
    def distance(self, p="const Vector&"):
        "Compute the minimum distance to a point."
        return "double"

    @PYB11const
    def closestPoint(self, p="const Vector&"):
        "Find the point in the polygon closest to the given point."
        return "Vector"

    @PYB11const
    def convex(self, tol = "double"):
        "Test if the polygon is convex"
        return "bool"

    def setBoundingBox(self):
        "Set the internal bounding box"
        return "void"

    @PYB11const
    def facetArea(self, facetID="const unsigned"):
        return "double"

    @PYB11const
    def facetAreaNormal(self, facetID="const unsigned"):
        return "Vector"

    @PYB11const
    def facetSubVolume(self, facetID="const unsigned"):
        "Decompose the polygon into triangles for each facet"
        return "Polygon"

    def transform(self, t="const Tensor&"):
        "Apply a general transformation tensor"
        return "Polygon&"

    #...........................................................................
    # Operators
    def __iadd__(self, rhs="Vector()"):
        return

    def __isub__(self, rhs="Vector()"):
        return

    def __add__(self, rhs="Vector()"):
        return

    def __sub__(self, rhs="Vector()"):
        return

    def __imul__(self, rhs="double()"):
        return
    
    def __itruediv__(self, rhs="double()"):
        return
    
    def __mul__(self, rhs="double()"):
        return
    
    def __truediv__(self, rhs="double()"):
        return
    
    def __eq__(self):
        return

    def __ne__(self):
        return

    #...........................................................................
    # Properties
    centroid = PYB11property("Vector")
    vertices = PYB11property("const std::vector<Vector>&", returnpolicy="reference_internal")
    facets = PYB11property("const std::vector<Facet>&", returnpolicy="reference_internal")
    facetVertices = PYB11property("std::vector<std::vector<unsigned> >",
                                  doc="Spit out a vector<vector<unsigned> > that encodes the facets.")
    vertexUnitNorms = PYB11property("const std::vector<Vector>&", returnpolicy="reference_internal")
    vertexFacetConnectivity = PYB11property("const std::vector<std::vector<unsigned> >&", returnpolicy="reference_internal")
    facetFacetConnectivity = PYB11property("const std::vector<std::vector<unsigned> >&", returnpolicy="reference_internal")
    xmin = PYB11property("const Vector&", returnpolicy="reference_internal")
    xmax = PYB11property("const Vector&", returnpolicy="reference_internal")
    volume = PYB11property("double")
