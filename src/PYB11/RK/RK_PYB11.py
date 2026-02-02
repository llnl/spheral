"""
Spheral RK module.

Provides reproducing kernel infrastructure.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from RKFieldNames import *
from RKCoefficients import *
from RKCorrections import *
from RKUtilities import *
from ReproducingKernelMethods import *
from ReproducingKernel import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"RK/RKCoefficients.hh"',
                  '"RK/RKCorrections.hh"',
                  '"RK/RKUtilities.hh"',
                  '"RK/RKFieldNames.hh"',
                  '"RK/ReproducingKernelMethods.hh"',
                  '"RK/ReproducingKernel.hh"',
                  '"RK/computeRKVolumes.hh"',
                  '"RK/computeOccupancyVolume.hh"',
                  '"RK/computeRKSumVolume.hh"',
                  '"RK/computeHullVolume.hh"',
                  '"RK/computeHullVolumes.hh"',
                  '"RK/computeHVolumes.hh"',
                  '"RK/computeOccupancyVolume.hh"',
                  '"RK/interpolateRK.hh"',
                  '"RK/gradientRK.hh"',
                  '"RK/hessianRK.hh"',
                  '"DataBase/State.hh"',
                  '"DataBase/StateDerivatives.hh"',
                  '"FileIO/FileIO.hh"',
                  '"Boundary/Boundary.hh"',
                  '"Utilities/NodeCoupling.hh"',
                  '<sstream>',
                  '<variant>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Preamble
#-------------------------------------------------------------------------------
PYB11preamble += """
namespace {
template<typename Dimension>
struct AppendFieldLists {
    
    typedef std::vector<std::variant<FieldList<Dimension, typename Dimension::Scalar>,
                                     FieldList<Dimension, typename Dimension::Vector>,
                                     FieldList<Dimension, typename Dimension::Tensor>,
                                     FieldList<Dimension, typename Dimension::SymTensor>,
                                     FieldList<Dimension, typename Dimension::ThirdRankTensor>>> FieldListArray;

    py::list& pylist;

    AppendFieldLists(py::list& pylist_):
        pylist(pylist_) {}

    template<typename FieldListType>
    inline
    void operator()(FieldListType& x) const {
        pylist.append(x);
    }
};
}
"""

#-------------------------------------------------------------------------------
# enums
#-------------------------------------------------------------------------------
RKOrder = PYB11enum(("ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder", "SexticOrder", "SepticOrder"),
                    export_values = True,
                    doc = "Selection of RK correction orders")
RKVolumeType = PYB11enum(("RKMassOverDensity", "RKSumVolume", "RKVoronoiVolume", "RKHullVolume", "HVolume"),
                         export_values = True,
                         doc = "Options for RK volume algorithms")

#-------------------------------------------------------------------------------
# Methods
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeRKVolumes(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                     W = "const TableKernel<%(Dimension)s>&",
                     position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                     mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                     massDensity = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                     H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                     damage = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                     facetedBoundaries = "const std::vector<typename %(Dimension)s::FacetedVolume>&",
                     facetedHoles = "const std::vector<std::vector<typename %(Dimension)s::FacetedVolume>>&",
                     boundaryConditions = "const std::vector<Boundary<%(Dimension)s>*>&",
                     volumeType = "const RKVolumeType",
                     surfacePoint = "FieldList<%(Dimension)s, int>&",
                     deltaCentroid = "FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                     etaVoidPoints = "FieldList<%(Dimension)s, std::vector<typename %(Dimension)s::Vector>>&",
                     cells = "FieldList<%(Dimension)s, typename %(Dimension)s::FacetedVolume>&",
                     cellFaceFlags = "FieldList<%(Dimension)s, std::vector<CellFaceFlag>>&",
                     volume = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute RK volumes"
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeRKSumVolume(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                       W = "const TableKernel<%(Dimension)s>&",
                       position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                       mass = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                       H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                       vol = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the RK mass density summation."
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeOccupancyVolume(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                           W = "const TableKernel<%(Dimension)s>&",
                           position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                           H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                           vol = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the occupancy volume per point"
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeHullVolume(position = "const FieldList<%(Dimension)s, %(Dimension)s::Vector>&",
                      H = "const FieldList<%(Dimension)s, %(Dimension)s::SymTensor>&",
                      connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                      clipToVoronoi = "const bool",
                      surfacePoint = "FieldList<%(Dimension)s, int>&",
                      vol = "FieldList<%(Dimension)s, %(Dimension)s::Scalar>&",
                      cells = "FieldList<%(Dimension)s, %(Dimension)s::FacetedVolume>&"):
    "Compute the volume per point based on convex hulls."
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeHullVolumes(connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                       kernelExtent = "const typename %(Dimension)s::Scalar",
                       position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                       H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                       volume = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute the volume per point based on convex hulls."
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
def computeHVolumes(nPerh = "const typename %(Dimension)s::Scalar",
                    H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                    volume = "FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&"):
    "Compute a volume per point based on the local H tensor."
    return "void"

#-------------------------------------------------------------------------------
@PYB11template("Dimension", "DataType")
@PYB11implementation("""[](const FieldList<%(Dimension)s, %(DataType)s>& fieldList,
                           const FieldList<%(Dimension)s, %(Dimension)s::Vector>& position,
                           const FieldList<%(Dimension)s, %(Dimension)s::Scalar>& weight,
                           const FieldList<%(Dimension)s, %(Dimension)s::SymTensor>& H,
                           const ConnectivityMap<%(Dimension)s>& connectivityMap,
                           const ReproducingKernel<%(Dimension)s>& WR,
                           const FieldList<%(Dimension)s, RKCoefficients<%(Dimension)s>>& corrections,
                           const NodeCoupling& nodeCoupling) {
                               std::vector<std::variant<FieldList<%(Dimension)s, %(Dimension)s::Scalar>,
                                                        FieldList<%(Dimension)s, %(Dimension)s::Vector>,
                                                        FieldList<%(Dimension)s, %(Dimension)s::Tensor>,
                                                        FieldList<%(Dimension)s, %(Dimension)s::SymTensor>,
                                                        FieldList<%(Dimension)s, %(Dimension)s::ThirdRankTensor>>> fieldLists;
                               fieldLists.emplace_back(fieldList);
                               auto flvec = interpolateRK(fieldLists,
                                                          position,
                                                          weight,
                                                          H,
                                                          connectivityMap,
                                                          WR,
                                                          corrections,
                                                          nodeCoupling);
                               CHECK(flvec.size() == 1);
                               FieldList<%(Dimension)s, %(DataType)s> result(std::get<FieldList<%(Dimension)s, %(DataType)s>>(flvec[0]));
                               result.copyFields();
                               return result;
                           }""")
@PYB11cppname("interpolateRK")
def interpolateRK1(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
                   position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                   weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                   H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                   connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                   WR = "const ReproducingKernel<%(Dimension)s>&",
                   corrections = "const FieldList<%(Dimension)s, RKCoefficients<%(Dimension)s>>&",
                   nodeCoupling = ("const NodeCoupling&", "NodeCoupling()")):
    "Compute the RK interpolation at each point for a single FieldList."
    return "FieldList<%(Dimension)s, %(DataType)s>"

#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11implementation("""[](const py::list& pyFieldLists,
                           const FieldList<%(Dimension)s, %(Dimension)s::Vector>& position,
                           const FieldList<%(Dimension)s, %(Dimension)s::Scalar>& weight,
                           const FieldList<%(Dimension)s, %(Dimension)s::SymTensor>& H,
                           const ConnectivityMap<%(Dimension)s>& connectivityMap,
                           const ReproducingKernel<%(Dimension)s>& WR,
                           const FieldList<%(Dimension)s, RKCoefficients<%(Dimension)s>>& corrections,
                           const NodeCoupling& nodeCoupling) {
                               std::vector<std::variant<FieldList<%(Dimension)s, %(Dimension)s::Scalar>,
                                                        FieldList<%(Dimension)s, %(Dimension)s::Vector>,
                                                        FieldList<%(Dimension)s, %(Dimension)s::Tensor>,
                                                        FieldList<%(Dimension)s, %(Dimension)s::SymTensor>,
                                                        FieldList<%(Dimension)s, %(Dimension)s::ThirdRankTensor>>> fieldLists;
                               for (auto& pyfl: pyFieldLists) {
                                   try {
                                           auto fl = pyfl.cast<FieldList<%(Dimension)s, %(Dimension)s::Scalar>&>();
                                           fieldLists.emplace_back(fl);
                                   } catch (py::cast_error& e) {
                                       try {
                                               auto fl = pyfl.cast<FieldList<%(Dimension)s, %(Dimension)s::Vector>&>();
                                               fieldLists.emplace_back(fl);
                                       } catch (py::cast_error& e) {
                                           try {
                                                   auto fl = pyfl.cast<FieldList<%(Dimension)s, %(Dimension)s::Tensor>&>();
                                                   fieldLists.emplace_back(fl);
                                           } catch (py::cast_error& e) {
                                               try {
                                                       auto fl = pyfl.cast<FieldList<%(Dimension)s, %(Dimension)s::SymTensor>&>();
                                                       fieldLists.emplace_back(fl);
                                               } catch (py::cast_error& e) {
                                                   try {
                                                           auto fl = pyfl.cast<FieldList<%(Dimension)s, %(Dimension)s::ThirdRankTensor>&>();
                                                           fieldLists.emplace_back(fl);
                                                   } catch (py::cast_error& e) {
                                                       throw py::type_error("interpolateCRKSPH expects a list of only (Scalar,Vector,Tensor,SymTensor,ThirdRankTensor)FieldList.");
                                                   }
                                               }
                                           }
                                       }
                                   }
                               }
                               auto cppresult =  interpolateRK(fieldLists,
                                                               position,
                                                               weight,
                                                               H,
                                                               connectivityMap,
                                                               WR,
                                                               corrections,
                                                               nodeCoupling);
                               py::list result;
                               for (auto fl: cppresult) std::visit(AppendFieldLists<%(Dimension)s>(result), fl);
                               return result;
                           }""")
def interpolateRK(fieldLists = "py::list&",
                  position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                  weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
                  H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
                  connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
                  WR = "const ReproducingKernel<%(Dimension)s>&",
                  corrections = "const FieldList<%(Dimension)s, RKCoefficients<%(Dimension)s>>&",
                  nodeCoupling = ("const NodeCoupling&", "NodeCoupling()")):
    "Compute the RK interpolation at each point for a list of FieldLists."
    return "py::list"

#-------------------------------------------------------------------------------
@PYB11template("Dimension", "DataType")
def gradientRK(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
               position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
               weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
               H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
               connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
               WR = "const ReproducingKernel<%(Dimension)s>&",
               corrections = "const FieldList<%(Dimension)s, RKCoefficients<%(Dimension)s>>&",
               nodeCoupling = ("const NodeCoupling&", "NodeCoupling()")):
    "Compute the RK gradient at each point for a FieldList."
    return "FieldList<%(Dimension)s, typename MathTraits<%(Dimension)s, %(DataType)s>::GradientType>"

#-------------------------------------------------------------------------------
@PYB11template("Dimension", "DataType")
@PYB11pycppname("gradientRK")
def gradientRK2(fieldList = "const FieldList<%(Dimension)s, std::vector<%(DataType)s>>&",
               position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
               weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
               H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
               connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
               WR = "const ReproducingKernel<%(Dimension)s>&",
               corrections = "const FieldList<%(Dimension)s, RKCoefficients<%(Dimension)s>>&",
               nodeCoupling = ("const NodeCoupling&", "NodeCoupling()")):
    "Compute the RK gradient at each point for a FieldList."
    return "FieldList<%(Dimension)s, std::vector<typename MathTraits<%(Dimension)s, %(DataType)s>::GradientType>>"

#-------------------------------------------------------------------------------
@PYB11template("Dimension", "DataType")
def hessianRK(fieldList = "const FieldList<%(Dimension)s, %(DataType)s>&",
              position = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
              weight = "const FieldList<%(Dimension)s, typename %(Dimension)s::Scalar>&",
              H = "const FieldList<%(Dimension)s, typename %(Dimension)s::SymTensor>&",
              connectivityMap = "const ConnectivityMap<%(Dimension)s>&",
              WR = "const ReproducingKernel<%(Dimension)s>&",
              corrections = "const FieldList<%(Dimension)s, RKCoefficients<%(Dimension)s>>&",
              nodeCoupling = ("const NodeCoupling&", "NodeCoupling()")):
    "Compute the RK hessian at each point for a FieldList."
    return "FieldList<%(Dimension)s, typename MathTraits<%(Dimension)s, %(DataType)s>::HessianType>"

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
RKCorrections%(ndim)id = PYB11TemplateClass(RKCorrections, template_parameters="%(Dimension)s")
ReproducingKernelMethods%(ndim)id = PYB11TemplateClass(ReproducingKernelMethods, template_parameters="%(Dimension)s")
ReproducingKernel%(ndim)id = PYB11TemplateClass(ReproducingKernel, template_parameters="%(Dimension)s")
RKCoefficients%(ndim)id = PYB11TemplateClass(RKCoefficients, template_parameters="%(Dimension)s")

interpolateRK%(ndim)id = PYB11TemplateFunction(interpolateRK, template_parameters="Dim<%(ndim)i>", pyname="interpolateRK")
computeRKVolumes%(ndim)id = PYB11TemplateFunction(computeRKVolumes, template_parameters="%(Dimension)s")
computeRKSumVolume%(ndim)id = PYB11TemplateFunction(computeRKSumVolume, template_parameters="%(Dimension)s")
computeOccupancyVolume%(ndim)id = PYB11TemplateFunction(computeOccupancyVolume, template_parameters="%(Dimension)s")
computeHullVolume%(ndim)id = PYB11TemplateFunction(computeHullVolume, template_parameters="%(Dimension)s")
computeHullVolumes%(ndim)id = PYB11TemplateFunction(computeHullVolumes, template_parameters="%(Dimension)s")
computeHVolumes%(ndim)id = PYB11TemplateFunction(computeHVolumes, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

    # RK interpolation
    for element in ("Dim<%i>::Scalar" % ndim,
                    "Dim<%i>::Vector" % ndim,
                    "Dim<%i>::Tensor" % ndim,
                    "Dim<%i>::SymTensor" % ndim,
                    "Dim<%i>::ThirdRankTensor" % ndim):
        exec('''
interpolateRK%(label)s = PYB11TemplateFunction(interpolateRK1, template_parameters=("%(Dimension)s", "%(element)s"), pyname="interpolateRK")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "element"   : element,
       "label"     : PYB11mangle(element)})

    # RK gradient
    for element in ("Dim<%i>::Scalar" % ndim,
                    "Dim<%i>::Vector" % ndim):
        exec('''
gradientRK%(label)s = PYB11TemplateFunction(gradientRK, template_parameters=("%(Dimension)s", "%(element)s"), pyname="gradientRK")
hessianRK%(label)s = PYB11TemplateFunction(hessianRK, template_parameters=("%(Dimension)s", "%(element)s"), pyname="hessianRK")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">",
       "element"   : element,
       "label"     : PYB11mangle(element)})

    # RKUtilities
    for num, correctionOrder in enumerate(("ZerothOrder", "LinearOrder", "QuadraticOrder", "CubicOrder", "QuarticOrder", "QuinticOrder", "SexticOrder", "SepticOrder")):
        exec(''' 
RKUtilities%(ndim)sd%(num)s = PYB11TemplateClass(RKUtilities, template_parameters=("%(Dimension)s", "RKOrder::%(correctionOrder)s"))
''' % {"ndim"            : ndim,
       "Dimension"       : "Dim<" + str(ndim) + ">",
       "num"             : num,
       "correctionOrder" : correctionOrder})

    # RK interpolation
