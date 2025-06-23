"""
Spheral Utilities module.

A bunch of utility methods used throughout Spheral.  Unfortunately this has become
a bit of a grab bag of math, geometry, infrastructure, and assorted functions with
no real relation.  Should probably revisit and categorize this stuff to other
modules more effectively.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Utilities/setGlobalFlags.hh"',
                  '"Utilities/packElement.hh"',
                  '"boost/math/special_functions/legendre.hpp"',
                  '"Utilities/BuildData.hh"',
                  '"Utilities/Functors.hh"',
                  '"Utilities/erff.hh"',
                  '"Utilities/newtonRaphson.hh"',
                  '"Utilities/bisectRoot.hh"',
                  '"Utilities/simpsonsIntegration.hh"',
                  '"Utilities/globalNodeIDs.hh"',
                  '"Utilities/rotationMatrix.hh"',
                  '"Utilities/iterateIdealH.hh"',
                  '"Utilities/nodeOrdering.hh"',
                  '"Utilities/mortonOrderIndices.hh"',
                  '"Utilities/peanoHilbertOrderIndices.hh"',
                  '"Utilities/boundingBox.hh"',
                  '"Utilities/globalBoundingVolumes.hh"',
                  '"Utilities/testBoxIntersection.hh"',
                  '"Utilities/lineSegmentIntersections.hh"',
                  '"Utilities/pointDistances.hh"',
                  '"Utilities/segmentIntersectEdges.hh"',
                  '"Utilities/integrateThroughMeshAlongSegment.hh"',
                  '"Utilities/numberDensity.hh"',
                  '"Utilities/planarReflectingOperator.hh"',
                  '"Utilities/KeyTraits.hh"',
                  '"Utilities/pointInPolygon.hh"',
                  '"Utilities/pointOnPolygon.hh"',
                  '"Utilities/pointInPolyhedron.hh"',
                  '"Utilities/pointOnPolyhedron.hh"',
                  '"Utilities/refinePolyhedron.hh"',
                  '"Utilities/overlayRemapFields.hh"',
                  '"Utilities/computeShepardsInterpolation.hh"',
                  '"Utilities/clipFacetedVolume.hh"',
                  '"Utilities/DomainNode.hh"',
                  '"Utilities/NodeCoupling.hh"',
                  '"Utilities/LinearInterpolator.hh"',
                  '"Utilities/QuadraticInterpolator.hh"',
                  '"Utilities/CubicHermiteInterpolator.hh"',
                  '"Utilities/XYInterpolator.hh"',
                  '"Utilities/BiLinearInterpolator.hh"',
                  '"Utilities/BiQuadraticInterpolator.hh"',
                  '"Utilities/BiCubicInterpolator.hh"',
                  '"Utilities/uniform_random.hh"',
                  '"Utilities/Timer.hh"',
                  '"Distributed/Communicator.hh"',
                  '"adiak.hpp"',
                  '<algorithm>']

#-------------------------------------------------------------------------------
# Preamble
#-------------------------------------------------------------------------------
PYB11preamble += """
namespace Spheral {

inline void spheral_adiak_init() {
  adiak::init((void*) Communicator::comm_ptr());
  // Always collect some curated default adiak information
  adiak::adiakversion();
  adiak::user();
  adiak::uid();
  adiak::launchdate();
  adiak::workdir();
  adiak::hostname();
  adiak::clustername();
  adiak::walltime();
  adiak::cputime();
  adiak::jobsize();
  adiak::numhosts();
  adiak::hostlist();
  adiak::mpi_library_version();
}

enum adiak_categories {
unset = 0,
all,
general,
performance,
control
};
}
"""

PYB11modulepreamble = """
TIME_PHASE_BEGIN("main");
Spheral::spheral_adiak_init();

// Call these routines when module is exited
auto atexit = py::module_::import("atexit");
atexit.attr("register")(py::cpp_function([]() {
   TIME_PHASE_END("main");
   adiak::fini();
   if (Spheral::TimerMgr::is_started()) {
      Spheral::TimerMgr::fini();
   } else {
      Communicator::finalize();
   }
}));
"""

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Standard functions
#-------------------------------------------------------------------------------
def setGlobalFlags():
    return

#-------------------------------------------------------------------------------
# Instantiate types and add dimension dependent functions.
#-------------------------------------------------------------------------------
from SpheralFunctor import *
from KeyTraits import *
from DomainNode import *
from NodeCoupling import *
from LinearInterpolator import *
from QuadraticInterpolator import *
from CubicHermiteInterpolator import *
from XYInterpolator import *
from BiLinearInterpolator import *
from BiQuadraticInterpolator import *
from BiCubicInterpolator import *
from uniform_random import *
from BuildData import *
from Adiak import *
from TimerMgr import *

ScalarScalarFunctor = PYB11TemplateClass(SpheralFunctor, template_parameters=("double", "double"))
ScalarPairScalarFunctor = PYB11TemplateClass(SpheralFunctor, template_parameters=("double", "std::pair<double,double>"))
ScalarScalarScalarFunctor = PYB11TemplateClass(Spheral2ArgFunctor, template_parameters=("double", "double", "double"))
ScalarScalarSymTensor2dFunctor = PYB11TemplateClass(Spheral2ArgFunctor, template_parameters=("double", "double", "Dim<2>::SymTensor"))

@PYB11template("Vector")
@PYB11implementation("[](std::vector<%(Vector)s>& positions) { %(Vector)s xmin, xmax; boundingBox(positions, xmin, xmax); return py::make_tuple(xmin, xmax); }")
def boundingBoxVec(positions = "std::vector<%(Vector)s>&"):
    "Minimum (axis-aligned) bounding box for a collection of %(Vector)s"
    return "py::tuple"

@PYB11template("Dimension")
@PYB11implementation("[](const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>& positions, const bool useGhosts) { %(Dimension)s::Vector xmin, xmax; boundingBox(positions, xmin, xmax, useGhosts); return py::make_tuple(xmin, xmax); }")
def boundingBoxFL(positions = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                  useGhosts = "const bool"):
    "Minimum (axis-aligned) bounding box for a FieldList<%(Dimension)s::Vector>"
    return "py::tuple"

@PYB11template("Dimension")
@PYB11implementation("[](const DataBase<%(Dimension)s>& dataBase) { %(Dimension)s::FacetedVolume nodeVolume, sampleVolume; globalBoundingVolumes(dataBase, nodeVolume, sampleVolume); return py::make_tuple(nodeVolume, sampleVolume); }")
def globalBoundingVolumes(dataBase = "const DataBase<%(Dimension)s>&"):
    "Return the bounding FacetedVolumes for the positions and (positions+h_extent) in the DataBase"
    return "py::tuple"

@PYB11template("Vector")
def collinear(a = "const %(Vector)s&",
              b = "const %(Vector)s&",
              c = "const %(Vector)s&",
              tol = ("const double", "1.0e-10")):
    "Test if points a, b, & c are collinear.   You should scale the tolerance handed in here!"
    return "bool"

@PYB11template("Vector")
def between(a = "const %(Vector)s&",
            b = "const %(Vector)s&",
            c = "const %(Vector)s&",
            tol = ("const double", "1.0e-10")):
    "Test if point c is between a & b.   You should scale the tolerance handed in here!"
    return "bool"

@PYB11template("Vector")
@PYB11cppname("segmentSegmentIntersection")
def segmentSegmentIntersectionTest(a0 = "const %(Vector)s&",
                                   a1 = "const %(Vector)s&",
                                   b0 = "const %(Vector)s&",
                                   b1 = "const %(Vector)s&",
                                   tol = ("const double", "1.0e-10")):
    """Check if two line segments intersect (does not compute intersection point).
The line segments are characterized by their endpoints: a_seg = (a0, a1)
                                                        b_seg = (b0, b1)
The advantage of this method is that it should be faster if you don't need
the actual intersection.
Based on stuff from "Computational Geometry in C", Joseph O'Rourke"""
    return "bool"

@PYB11template("Dimension")
@PYB11cppname("numGlobalNodes")
def numGlobalNodesNL(nodes = "const NodeList<%(Dimension)s>&"):
    "Total number of nodes in the NodeList across all domains"
    return "int"

@PYB11template("Dimension")
@PYB11cppname("numGlobalNodes")
def numGlobalNodesDB(dataBase = "const DataBase<%(Dimension)s>&"):
    "Total number of nodes in the DataBase across all domains"
    return "int"

@PYB11template("Dimension")
@PYB11cppname("globalNodeIDs")
def globalNodeIDsNL(dataBase = "const NodeList<%(Dimension)s>&"):
    """Compute a unique set of global node IDs for the given NodeList, and return
the set of them on this process."""
    return "Field<%(Dimension)s, int>"

@PYB11template("Dimension")
@PYB11cppname("globalNodeIDs")
def globalNodeIDsDB(dataBase = "const DataBase<%(Dimension)s>&"):
    """Compute a unique set of global node IDs for all nodes across all NodeLists in
a DataBase, returning the result as a FieldList<int>."""
    return "FieldList<%(Dimension)s, int>"

@PYB11template("Dimension")
def iterateIdealH(dataBase = "DataBase<%(Dimension)s>&",
                  packages = "std::vector<Physics<%(Dimension)s>*>&",
                  boundaries = "const std::vector<Boundary<%(Dimension)s>*>&",
                  maxIterations = ("const int", "100"),
                  tolerance = ("const double", "1.0e-10"),
                  nPerhForIteration = ("const double", "0.0"),
                  sphericalStart = ("const bool", "false"),
                  fixDeterminant = ("const bool", "false")):
    """Iterate the ideal H algorithm to converge on a new H field.
This routine replaces the H field in place."""
    return "void"

@PYB11template("Dimension", "DataType")
def nodeOrdering(criteria = "const FieldList<%(Dimension)s, %(DataType)s>&"):
    """Compute the order that a given set of nodes should be stepped through
given a FieldList of things to sort them by.
The FieldList returned is the one to N indexing corresponding to sorting the 
input in increasing order."""
    return "FieldList<%(Dimension)s, int>"

@PYB11template("Dimension")
@PYB11cppname("mortonOrderIndices")
def mortonOrderIndices_pos(positions = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&"):
    "Return the Morton (Z) ordering for the positions"
    return "FieldList<%(Dimension)s, typename KeyTraits::Key>"

@PYB11template("Dimension")
@PYB11cppname("mortonOrderIndices")
def mortonOrderIndices_db(dataBase = "const DataBase<%(Dimension)s>&"):
    "Return the Morton (Z) ordering for the positions in the DataBase"
    return "FieldList<%(Dimension)s, typename KeyTraits::Key>"

@PYB11template("Dimension")
@PYB11cppname("mortonOrderIndices")
def mortonOrderIndices_mask(dataBase = "const DataBase<%(Dimension)s>&",
                            mask = "const FieldList<%(Dimension)s, int>&"):
    """Return the Morton (Z) ordering for the positions in the DataBase
This version allows some nodes to be screened with a mask"""
    return "FieldList<%(Dimension)s, typename KeyTraits::Key>"

@PYB11template("Dimension")
@PYB11cppname("peanoHilbertOrderIndices")
def peanoHilbertOrderIndices_pos(positions = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&"):
    "Return the Peano-Hilbert ordering for the positions"
    return "FieldList<%(Dimension)s, typename KeyTraits::Key>"

@PYB11template("Dimension")
@PYB11cppname("peanoHilbertOrderIndices")
def peanoHilbertOrderIndices_db(dataBase = "const DataBase<%(Dimension)s>&"):
    "Return the Peano-Hilbert ordering for the positions in the DataBase"
    return "FieldList<%(Dimension)s, typename KeyTraits::Key>"

@PYB11template("Dimension")
def numberDensity(dataBase = "const DataBase<%(Dimension)s>&",
                  W = "const TableKernel<%(Dimension)s>&"):
    """Return the number density for all nodes in the DataBase.
Note this method assumes that all boundary conditions, ghost nodes, and 
neighbor information is up to date."""
    return "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>"

@PYB11template("Dimension", "Value")
def integrateThroughMeshAlongSegment(values = "const std::vector<std::vector<%(Value)s>>&",
                                     xmin = "const typename %(Dimension)s::Vector&",
                                     xmax = "const typename %(Dimension)s::Vector&",
                                     ncells = "const std::vector<unsigned>&",
                                     s0 = "const typename %(Dimension)s::Vector&",
                                     s1 = "const typename %(Dimension)s::Vector&"):
    """Return the result of integrating a quantity along a line segment.
The quantity here is assumed to be represented a values in a vector<%(Value)s>,
where the vector<%(Value)s> is the value of the quantity in a series of cartesian
cells whose box is defined by by xmin, xmax, and ncells.

We actually pass in a vector<vector<%(Value)s>>, which is a progressively refined
(by factors of 2 in each dimesion) representation of the data.  The idea is that
we use the finest level with a non-zero value for the value."""
    return "%(Value)s"

@PYB11template("Dimension", "DataType")
def computeShepardsInterpolation(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
                                 connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                                 W = "const TableKernel<%(Dimension)s>&",
                                 position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                                 H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                                 weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the Shepard interpolation of a FieldList."
    return "FieldList<%(Dimension)s, %(DataType)s>"

#...............................................................................
# Instantiate stuff for the dimensions Spheral is building
for ndim in dims:
    exec('''
# Functors
VectorScalarFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Vector)s", "double"))
VectorVectorFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Vector)s", "%(Vector)s"))
VectorPairScalarFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Vector)s", "std::pair<double,double>"))
SizetSizetSymTensorSymTensorSymTensorFunctor%(ndim)id = PYB11TemplateClass(Spheral4ArgFunctor, template_parameters=("size_t", "size_t", "%(SymTensor)s", "%(SymTensor)s", "%(SymTensor)s"))

# boundingVolumes
boundingBoxVec%(ndim)id = PYB11TemplateFunction(boundingBoxVec, template_parameters="%(Vector)s", pyname="boundingBox")
boundingBoxFL%(ndim)id = PYB11TemplateFunction(boundingBoxFL, template_parameters="%(Dimension)s", pyname="boundingBox")
globalBoundingVolumes%(ndim)id = PYB11TemplateFunction(globalBoundingVolumes, template_parameters="%(Dimension)s", pyname="globalBoundingVolumes")

collinear%(ndim)id = PYB11TemplateFunction(collinear, template_parameters="%(Vector)s", pyname="collinear")
between%(ndim)id = PYB11TemplateFunction(between, template_parameters="%(Vector)s", pyname="between")
segmentSegmentIntersectionTest%(ndim)id = PYB11TemplateFunction(segmentSegmentIntersectionTest, template_parameters="%(Vector)s", pyname="segmentSegmentIntersectionTest")

DomainNode%(ndim)id = PYB11TemplateClass(DomainNode, template_parameters="%(Dimension)s")
vector_of_DomainNode%(ndim)id = PYB11_bind_vector("DomainNode<%(Dimension)s>", opaque=True, local=False)

#...............................................................................
# segment-segment intersections (return intersect)
@PYB11pycppname("segmentSegmentIntersection")
@PYB11implementation("""[](const %(Vector)s& a0,
                           const %(Vector)s& a1,
                           const %(Vector)s& b0,
                           const %(Vector)s& b1,
                           const double tol) { %(Vector)s result1, result2; auto flag = segmentSegmentIntersection(a0, a1, b0, b1, result1, result2, tol); return py::make_tuple(flag, result1, result2); }""")
def segmentSegmentIntersection%(ndim)id(a0 = "const %(Vector)s&",
                                        a1 = "const %(Vector)s&",
                                        b0 = "const %(Vector)s&",
                                        b1 = "const %(Vector)s&",
                                        tol = ("const double", "1.0e-10")):
    """Intersection of two line segments.
The line segments are characterized by their endpoints: a_seg = (a0, a1)
                                                        b_seg = (b0, b1)
Return values are a tuple(char, Vector, Vector)
       The char is a code characterizing the intersection:
             "e" -> The segments colinearly overlap (edge)
             "v" -> The endpoint of one segment lies on the other
             "1" -> The segments intersect properly
             "0" -> The segments do not intersect
       The Vectors are the intersection points (if any)"""
    return "py::tuple"

#...............................................................................
# segment-segment distance
@PYB11pycppname("segmentSegmentDistance")
def segmentSegmentDistance%(ndim)id(a0 = "const %(Vector)s&",
                                    a1 = "const %(Vector)s&",
                                    b0 = "const %(Vector)s&",
                                    b1 = "const %(Vector)s&"):
    "Compute the minimum distance between two line segments. (a0,a1) -> (b0,b1)"
    return "double"

#...............................................................................
numGlobalNodes%(ndim)id = PYB11TemplateFunction(numGlobalNodesNL, template_parameters="%(Dimension)s")
numGlobalNodesAll%(ndim)id = PYB11TemplateFunction(numGlobalNodesDB, template_parameters="%(Dimension)s")
globalNodeIDs%(ndim)id = PYB11TemplateFunction(globalNodeIDsNL, template_parameters="%(Dimension)s")
globalNodeIDsAll%(ndim)id = PYB11TemplateFunction(globalNodeIDsDB, template_parameters="%(Dimension)s")

iterateIdealH%(ndim)id = PYB11TemplateFunction(iterateIdealH, template_parameters="%(Dimension)s")

nodeOrdering%(ndim)id = PYB11TemplateFunction(nodeOrdering, template_parameters=("%(Dimension)s", "KeyTraits::Key"))
mortonOrderIndices_pos%(ndim)id = PYB11TemplateFunction(mortonOrderIndices_pos, template_parameters="%(Dimension)s", pyname="mortonOrderIndices%(ndim)id")
mortonOrderIndices_db%(ndim)id = PYB11TemplateFunction(mortonOrderIndices_db, template_parameters="%(Dimension)s", pyname="mortonOrderIndices%(ndim)id")
mortonOrderIndices_mask%(ndim)id = PYB11TemplateFunction(mortonOrderIndices_mask, template_parameters="%(Dimension)s", pyname="mortonOrderIndices%(ndim)id")
peanoHilbertOrderIndices_pos%(ndim)id = PYB11TemplateFunction(peanoHilbertOrderIndices_pos, template_parameters="%(Dimension)s", pyname="peanoHilbertOrderIndices%(ndim)id")
peanoHilbertOrderIndices_db%(ndim)id = PYB11TemplateFunction(peanoHilbertOrderIndices_db, template_parameters="%(Dimension)s", pyname="peanoHilbertOrderIndices%(ndim)id")

numberDensity%(ndim)id = PYB11TemplateFunction(numberDensity, template_parameters="%(Dimension)s")
integrateThroughMeshAlongSegment%(ndim)id = PYB11TemplateFunction(integrateThroughMeshAlongSegment, template_parameters=("%(Dimension)s", "double"))

computeShepardsInterpolationScalar%(ndim)id = PYB11TemplateFunction(computeShepardsInterpolation, template_parameters=("%(Dimension)s", "double"))
computeShepardsInterpolationVector%(ndim)id = PYB11TemplateFunction(computeShepardsInterpolation, template_parameters=("%(Dimension)s", "%(Vector)s"))
computeShepardsInterpolationTensor%(ndim)id = PYB11TemplateFunction(computeShepardsInterpolation, template_parameters=("%(Dimension)s", "%(Tensor)s"))
computeShepardsInterpolationSymTensor%(ndim)id = PYB11TemplateFunction(computeShepardsInterpolation, template_parameters=("%(Dimension)s", "%(SymTensor)s"))
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "Vector"    : "Dim<" + str(ndim) + ">::Vector",
       "Tensor"    : "Dim<" + str(ndim) + ">::Tensor",
       "SymTensor" : "Dim<" + str(ndim) + ">::SymTensor"})

#-------------------------------------------------------------------------------
# Some methods are always instantiated
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def planarReflectingOperator(nhat = "const %(Dimension)s::Vector&"):
    "Generate the reflection operator for the (outward facing) normal"
    return "%(Dimension)s::Tensor"

@PYB11template("Dimension")
@PYB11cppname("planarReflectingOperator")
def planarReflectingOperator1(plane = "const GeomPlane<%(Dimension)s>&"):
    "Generate the reflection operator for the given plane."
    return "%(Dimension)s::Tensor"

for ndim in range(1, 4):
    exec('''
@PYB11cppname("rotationMatrix")
def rotationMatrix%(ndim)id(runit = "const %(Vector)s&"):
    "Determine the rotation matrix to rotate into the frame such that the given vector is aligned with the x axis."
    return "%(Tensor)s"

@PYB11cppname("testBoxIntersection")
def testBoxIntersection%(ndim)id(xmin1 = "const %(Vector)s&",
                                 xmax1 = "const %(Vector)s&",
                                 xmin2 = "const %(Vector)s&",
                                 xmax2 = "const %(Vector)s&",
                                 tol = ("const double", "1.0e-10")):
    "Test if two boxes (specified by their min and max corners) intersect"
    return "bool"

@PYB11cppname("testPointInBox")
def testPointInBox%(ndim)id(point = "const %(Vector)s&",
                            xmin = "const %(Vector)s&",
                            xmax = "const %(Vector)s&",
                            tol = ("const double", "1.0e-10")):
    "Test if a point is contained in a box."
    return "bool"

planarReflectingOperator%(ndim)id = PYB11TemplateFunction(planarReflectingOperator, template_parameters="%(Dimension)s", pyname="planarReflectingOperator")
planarReflectingOperator1%(ndim)id = PYB11TemplateFunction(planarReflectingOperator1, template_parameters="%(Dimension)s", pyname="planarReflectingOperator")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "Vector"    : "Dim<" + str(ndim) + ">::Vector",
       "Tensor"    : "Dim<" + str(ndim) + ">::Tensor",
       "SymTensor" : "Dim<" + str(ndim) + ">::SymTensor",
       "Plane"     : "GeomPlane<Dim<" + str(ndim) + ">>"})

#-------------------------------------------------------------------------------
# Module functions
#-------------------------------------------------------------------------------
def erff(x = "double"):
    "You know, the error function."
    return "double"

@PYB11namespace("boost::math")
def legendre_p(l = "int",
               m = "int",
               x = "double"):
    "Compute the associated Legendre polynomial P^m_l(x)"
    return "double"

@PYB11cppname("bisectRoot<const PythonBoundFunctors::SpheralFunctor<double, double>>")
def bisectRoot(function = "const PythonBoundFunctors::SpheralFunctor<double, double>&",
               xmin = "double",
               xmax = "double",
               xaccuracy = ("double", "1.0e-15"),
               yaccuracy = ("double", "1.0e-10"),
               maxIterations = ("unsigned", "100u"),
               verbose = ("bool", "false")):
    """Bisection root finder.
Finds a root of 'function' in the range (x1, x2)"""
    return "double"

@PYB11cppname("newtonRaphson<const PythonBoundFunctors::SpheralFunctor<double, std::pair<double, double>>>")
def newtonRaphsonFindRoot(function = "const PythonBoundFunctors::SpheralFunctor<double, std::pair<double, double>>&",
                          x1 = "double",
                          x2 = "double",
                          xaccuracy = ("double", "1.0e-15"),
                          yaccuracy = ("double", "1.0e-15"),
                          maxIterations = ("int", "100")):
    """Newton-Raphson root finder.
Finds a root of 'function' in the range (x1, x2)"""
    return "double"

@PYB11cppname("simpsonsIntegration<const PythonBoundFunctors::SpheralFunctor<double, double>, double, double>")
def simpsonsIntegrationDouble(function = "const PythonBoundFunctors::SpheralFunctor<double, double>&",
                              x0 = "double",
                              x1 = "double",
                              numBins = "unsigned"):
    "Numerically integrate 'function' in the range (x0, x1) via Simpsons rule"
    return "double"

#-------------------------------------------------------------------------------
# packElement/unpackElement
#-------------------------------------------------------------------------------
@PYB11template("T")
@PYB11implementation("[](const %(T)s& x) -> py::bytes { std::vector<char> buf; packElement(x, buf); std::string strbuf; strbuf.resize(buf.size()); std::copy(buf.begin(), buf.end(), strbuf.begin()); return py::bytes(strbuf); }")
def packElement(x = "const %(T)s&"):
    "Serialize a %(T)s into a buffer of py::bytes"
    return "py::bytes"

@PYB11template("T")
@PYB11implementation("[](const std::string& x) -> %(T)s { %(T)s result; std::vector<char> buf(x.size()); std::copy(x.begin(), x.begin() + x.size(), buf.begin()); std::vector<char>::const_iterator itr = buf.begin(); std::vector<char>::const_iterator endItr = buf.end(); unpackElement(result, itr, endItr); return result; }")
def unpackElement(x = "const std::string&"):
    "Deserialize a %(T)s from a string buffer"
    return "%(T)s"

for type, suffix in (("int", "Int"),
                     ("unsigned", "Unsigned"),
                     ("double", "Double"),
                     ("uint32_t", "UL"),
                     ("uint64_t", "ULL"),
                     ("Dim<1>::FacetedVolume", "Box1d"),
                     ("Dim<2>::FacetedVolume", "Polygon"),
                     ("Dim<3>::FacetedVolume", "Polyhedron"),
                     ("std::vector<int>", "VectorOfInt"),
                     ("std::vector<unsigned>", "VectorOfUnsigned"),
                     ("std::vector<float>", "VectorOfFloat"),
                     ("std::vector<double>", "VectorOfDouble"),
                     ("std::vector<uint64_t>", "VectorOfULL"),
                     ("std::vector<std::vector<unsigned>>", "VectorOfVectorOfUnsigned")):
    exec('''
packElement%(suffix)s       = PYB11TemplateFunction(packElement,   template_parameters="%(type)s")   # Backwards compatible name
packElement2%(suffix)s      = PYB11TemplateFunction(packElement,   template_parameters="%(type)s", pyname="packElement")
unpackElement%(suffix)s     = PYB11TemplateFunction(unpackElement, template_parameters="%(type)s")
''' % {"type"   : type,
       "suffix" : suffix})

#-------------------------------------------------------------------------------
# Geometry operations only available in 2D/3D.
#-------------------------------------------------------------------------------
@PYB11template("Vector")
def closestPointOnSegment(p = "const %(Vector)s&",
                          a0 = "const %(Vector)s&",
                          a1 = "const %(Vector)s&"):
    "Find the point on a line segment (a0,a1) closest to point (p)."
    return "%(Vector)s"

@PYB11template("Vector")
@PYB11cppname("closestPointOnSegment")
def closestPointOnSegment1(p = "const %(Vector)s&",
                           a0 = "const %(Vector)s&",
                           a1 = "const %(Vector)s&",
                           result = "%(Vector)s&"):
    """Find the point on a line segment (a0,a1) closest to point (p).
This version return True if the closest point on the line is bounded by the segment, and False otherwise."""
    return "bool"

@PYB11template("Vector")
def pointPlaneDistance(point = "const %(Vector)s&",
                       origin = "const %(Vector)s&",
                       unitNormal = "const %(Vector)s&"):
    "Compute the distance between a point and a plane"
    return "double"

@PYB11template("Dimension")
def overlayRemapFields(boundaries = "const std::vector<Boundary<%(Dimension)s>*>&",
                       scalarDonorFields = "const std::vector<Field<%(Dimension)s, typename %(Dimension)s::Scalar>*>&",
                       vectorDonorFields = "const std::vector<Field<%(Dimension)s, typename %(Dimension)s::Vector>*>&",
                       tensorDonorFields = "const std::vector<Field<%(Dimension)s, typename %(Dimension)s::Tensor>*>&",
                       symTensorDonorFields = "const std::vector<Field<%(Dimension)s, typename %(Dimension)s::SymTensor>*>&",
                       scalarAcceptorFields = "std::vector<Field<%(Dimension)s, typename %(Dimension)s::Scalar>*>&",
                       vectorAcceptorFields = "std::vector<Field<%(Dimension)s, typename %(Dimension)s::Vector>*>&",
                       tensorAcceptorFields = "std::vector<Field<%(Dimension)s, typename %(Dimension)s::Tensor>*>&",
                       symTensorAcceptorFields = "std::vector<Field<%(Dimension)s, typename %(Dimension)s::SymTensor>*>&"):
    """Use geometric clipping to remap a set of conserved fields.
Currently only works single NodeList -> single NodeList, no boundaries."""
    return "void"

for ndim in (x for x in dims if x in (2, 3)):
    exec('''
closestPointOnSegment%(ndim)id = PYB11TemplateFunction(closestPointOnSegment, template_parameters="%(Vector)s", pyname="closestPointOnSegment")
closestPointOnSegment1%(ndim)id = PYB11TemplateFunction(closestPointOnSegment1, template_parameters="%(Vector)s", pyname="closestPointOnSegment")
pointPlaneDistance%(ndim)id = PYB11TemplateFunction(pointPlaneDistance, template_parameters="%(Vector)s", pyname="pointPlaneDistance")
overlayRemapFields%(ndim)id = PYB11TemplateFunction(overlayRemapFields, template_parameters="Dim<%(ndim)i>", pyname="overlayRemapFields")
''' % {"ndim"   : ndim,
       "Vector" : "Dim<" + str(ndim) + ">::Vector"})

#-------------------------------------------------------------------------------
# Geometry operations only available in 3D.
#-------------------------------------------------------------------------------
# segment-planar section intersection
@PYB11implementation("""[](const Dim<3>::Vector& s0,
                           const Dim<3>::Vector& s1,
                           const std::vector<Dim<3>::Vector>& pverts,
                           const double tol) {
                               Dim<3>::Vector intersect;
                               auto flag = segmentPlanarSectionIntersection(s0, s1, pverts, intersect, tol);
                               return py::make_tuple(flag, intersect);
                           }""")
def segmentPlanarSectionIntersection(s0 = "const Dim<3>::Vector&",
                                     s1 = "const Dim<3>::Vector&",
                                     pverts = "const std::vector<Dim<3>::Vector>&",
                                     tol = ("const double", "1.0e-8")):
    """Intersection of a line segment with a polygonal section of a plane.
The line segment is characterized by it's endpoints:     seg = (s0, s1)
The polygonal section of the plane is specified by
a series of points:                                    plane = pverts[0], pverts[1], ...
Note there must be at least 3 non-collinear points, and they
must be passed in in order to draw the polygon.

Return values are a tuple(char, Vector)
       The Vector is the intersection point (if any)
       The char is a code characterizing the intersection:
             "p" -> The segment lies in the plane (plane)
             "d" -> The pverts points do not define a unique plane (degenerate)
             "1" -> The segment intersects the plane properly
             "0" -> The segment does not intersect the plane"""
    return "py::tuple"

#...............................................................................
# segment-plane intersection (point-normal plane)
@PYB11implementation("""[](const Dim<3>::Vector& s0,
                           const Dim<3>::Vector& s1,
                           const Dim<3>::Vector& point,
                           const Dim<3>::Vector& normal,
                           const double tol) {
                               Dim<3>::Vector intersect;
                               auto flag = segmentPlaneIntersection(s0, s1, point, normal, intersect, tol);
                               return py::make_tuple(flag, intersect);
                           }""")
def segmentPlaneIntersection(s0 = "const Dim<3>::Vector&",
                             s1 = "const Dim<3>::Vector&",
                             point = "const Dim<3>::Vector&",
                             normal = "const Dim<3>::Vector&",
                             tol = ("const double", "1.0e-8")):
    """Intersection of a line segment with a plane.
The line segment is characterized by it's endpoints:     seg = (s0, s1)
The plane is characterized by a point in the plane and a unit normal: plane (point, normal)

Return values are a tuple<char, Vector>
       The Vector is the intersection point (if any)
       The char is a code characterizing the intersection:
             "p" -> The segment lies in the plane (plane)
             "d" -> The p points do not define a unique plane (degenerate)
             "1" -> The segment intersects the plane properly
             "0" -> The segment does not intersect the plane"""
    return "py::tuple"

#...............................................................................
# segment-plane intersection (three point plane)
@PYB11implementation("""[](const Dim<3>::Vector& s0,
                           const Dim<3>::Vector& s1,
                           const Dim<3>::Vector& p0,
                           const Dim<3>::Vector& p1,
                           const Dim<3>::Vector& p2,
                           const double tol) {
                               Dim<3>::Vector intersect;
                               auto flag = segmentPlaneIntersection(s0, s1, p0, p1, p2, intersect, tol);
                               return py::make_tuple(flag, intersect);
                           }""")
@PYB11pycppname("segmentPlaneIntersection")
def segmentPlaneIntersection1(s0 = "const Dim<3>::Vector&",
                              s1 = "const Dim<3>::Vector&",
                              p0 = "const Dim<3>::Vector&",
                              p1 = "const Dim<3>::Vector&",
                              p2 = "const Dim<3>::Vector&",
                              tol = ("const double", "1.0e-8")):
    """Intersection of a line segment with a plane.
The line segment is characterized by it's endpoints:     seg = (s0, s1)
The plane is characterized by three points in the plane: plane = (p0, p1, p2)

Return values are a tuple<char, Vector>
       The Vector is the intersection point (if any)
       The char is a code characterizing the intersection:
             "p" -> The segment lies in the plane (plane)
             "d" -> The p points do not define a unique plane (degenerate)
             "1" -> The segment intersects the plane properly
             "0" -> The segment does not intersect the plane"""
    return "py::tuple"

#...............................................................................
# segment-plane intersection True/False test
@PYB11pycppname("segmentPlaneIntersection")
def segmentPlaneIntersection2(a0 = "const Dim<3>::Vector&",
                              a1 = "const Dim<3>::Vector&",
                              vertices = "const std::vector<Dim<3>::Vector>&",
                              normal = "const Dim<3>::Vector&",
                              tol = ("const double", "1.0e-10")):
    """Check if a line segment intersects a polygonal planar section in 3D.
The line segment is characterized by its endpoints: a_seg = (a0, a1)
The planar section is characterized the a set of coplanar vertices: vertices.

The advantage of this method is that it should be faster if you don't need
the actual intersection.

This version takes the plane representation as a set of coplanar vertices."""
    return "bool"

#...............................................................................
# segment-plane intersection True/False test
@PYB11pycppname("segmentPlaneIntersection")
def segmentPlaneIntersection3(a0 = "const Dim<3>::Vector&",
                              a1 = "const Dim<3>::Vector&",
                              vertices = "const std::vector<Dim<3>::Vector>&",
                              ipoints = "const std::vector<unsigned>&",
                              normal = "const Dim<3>::Vector&",
                              tol = ("const double", "1.0e-10")):
    """Check if a line segment intersects a polygonal planar section in 3D.
The line segment is characterized by its endpoints: a_seg = (a0, a1)
The planar section is characterized the a set of coplanar vertices: vertices.

The advantage of this method is that it should be faster if you don't need
the actual intersection.

This version takes the plane representation as a subset of the vertices 
given by ipoints"""
    return "bool"

#...............................................................................
def closestPointOnPlane(p = "const Dim<3>::Vector&",
                        origin = "const Dim<3>::Vector&",
                        unitNormal = "const Dim<3>::Vector&"):
    "Find the point on a plane closest to the given point."
    return "Dim<3>::Vector"

#-------------------------------------------------------------------------------
# Polygon/Polyhedron containment
#-------------------------------------------------------------------------------
for volume, name, Vector in (("Dim<2>::FacetedVolume", "Polygon",    "Dim<2>::Vector"),
                             ("Dim<3>::FacetedVolume", "Polyhedron", "Dim<3>::Vector")):
    exec('''
def pointOn%(name)s(p = "const %(Vector)s&",
                    poly = "const %(volume)s&",
                    tol = ("const double", "1.0e-10")):
    "Test if a point is on the surface of a %(name)s"
    return "bool"

def pointIn%(name)s(p = "const %(Vector)s&",
                    poly = "const %(volume)s&",
                    countBoundary = ("const bool", "false"),
                    tol = ("const double", "1.0e-10")):
    "Test if a point is contained in a %(name)s"
    return "bool"

@PYB11pycppname("segmentIntersectEdges")
def segmentIntersectEdges%(name)s(a0 = "const %(Vector)s&",
                                  a1 = "const %(Vector)s&",
                                  poly = "const %(volume)s&",
                                  tol = ("const double", "1.0e-8")):
    "Test if a line segment (a0,a1) intersects any edges of a %(name)s"
    return "bool"
''' % {"volume" : volume,
       "name"   : name,
       "Vector" : Vector})

#...............................................................................
@PYB11pycppname("pointInPolygon")
def pointInPolygon1(p = "const Dim<2>::Vector&",
                    vertices = "const std::vector<Dim<2>::Vector>&"):
    "Test if point (p) is in polygon defined by (vertices)"
    return "bool"

#...............................................................................
@PYB11pycppname("pointInPolygon")
def pointInPolygon2(p = "const Dim<3>::Vector&",
                    vertices = "const std::vector<Dim<3>::Vector>&",
                    normal = "const Dim<3>::Vector&"):
    """3D point in polygon test -- the point and all vertices of the polygon must be coplanar.
  p        : point we're testing
  vertices : the vertices of the polygon
  normal   : normal to the plane of the polygon"""
    return "bool"

#...............................................................................
@PYB11pycppname("pointInPolygon")
def pointInPolygon3(p = "const Dim<3>::Vector&",
                    vertices = "const std::vector<Dim<3>::Vector>&",
                    ipoints = "const std::vector<unsigned>&",
                    normal = "const Dim<3>::Vector&",
                    countBoundary = ("const bool", "false"),
                    tol = ("const double", "1.0e-10")):
    """3D point in polygon test -- the point and all vertices of the polygon must be coplanar.
  p             : point we're testing
  vertices      : contains the vertices the polygon as a subset
  ipoints       : the subset of vertices that define the polygon
  normal        : normal to the plane of the polygon
  countBoundary : should points on the boundary count as internal?
  tol           : tolerance for on boundary"""
    return "bool"

#...............................................................................
def refinePolyhedron(poly0 = "const Dim<3>::FacetedVolume&",
                     numLevels = "int"):
    "Return a new Polyhedron based on refining an existing one a given number of levels."
    return "Dim<3>::FacetedVolume"

#-------------------------------------------------------------------------------
# PolyClipper utilities
#-------------------------------------------------------------------------------
@PYB11pycppname("clipFacetedVolume")
def clipFacetedVolume1(poly = "const Dim<2>::FacetedVolume&",
                       planes = "const std::vector<GeomPlane<Dim<2>>>&"):
    "Clip a polygon by a series of planes using R3D"
    return "Dim<2>::FacetedVolume"

@PYB11pycppname("clipFacetedVolume")
def clipFacetedVolume2(poly = "const Dim<3>::FacetedVolume&",
                       planes = "const std::vector<GeomPlane<Dim<3>>>&"):
    "Clip a polyhedron by a series of planes using R3D"
    return "Dim<3>::FacetedVolume"

#...............................................................................
@PYB11pycppname("clippedVolume")
def clippedVolume(poly = "const Dim<2>::FacetedVolume&",
                  planes = "const std::vector<GeomPlane<Dim<2>>>&"):
    "Return the volume of the clipped region."
    return "double"

@PYB11pycppname("clippedVolume")
def clippedVolume(poly = "const Dim<3>::FacetedVolume&",
                  planes = "const std::vector<GeomPlane<Dim<3>>>&"):
    "Return the volume of the clipped region."
    return "double"

#...............................................................................
for (value, label) in (("int", "Int"),
                       ("unsigned", "Unsigned"),
                       ("long", "Long"),
                       ("double", "Scalar"),
                       ("std::string", "String")):
    exec(f"""
adiak_value{label} = PYB11TemplateFunction(adiak_value, "{value}", pyname="adiak_value")
adiak_value2{label} = PYB11TemplateFunction(adiak_value2, "{value}", pyname="adiak_value")
""")
