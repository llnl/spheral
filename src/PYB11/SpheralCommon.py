#-------------------------------------------------------------------------------
# Common PYB11 initialization code for all Spheral modules.
#-------------------------------------------------------------------------------
from PYB11Generator import *

PYB11includes = ['"Geometry/Dimension.hh"',
                 '"Geometry/GeomPlane.hh"',
                 '"RK/RKCoefficients.hh"',
                 '"Field/Field.hh"',
                 '"Field/FieldList.hh"',
                 '"DataBase/DataBase.hh"',
                 '"DataOutput/registerWithRestart.hh"',
                 "<vector>",
                 "<map>",
                 "<set>",
                 "<string>",
                 '"polyclipper2d.hh"',
                 '"polyclipper3d.hh"']

PYB11preamble = """
typedef Spheral::GeomPlane<Spheral::Dim<1>> Plane1d;
typedef Spheral::Dim<1>::Vector Vector1d;
typedef Spheral::Dim<1>::Tensor Tensor1d;
typedef Spheral::Dim<1>::SymTensor SymTensor1d;
typedef Spheral::Dim<1>::ThirdRankTensor ThirdRankTensor1d;
typedef Spheral::Dim<1>::FourthRankTensor FourthRankTensor1d;
typedef Spheral::Dim<1>::FifthRankTensor FifthRankTensor1d;
typedef Spheral::Dim<1>::FacetedVolume FacetedVolume1d;
typedef Spheral::GeomPlane<Spheral::Dim<2>> Plane2d;
typedef Spheral::Dim<2>::Vector Vector2d;
typedef Spheral::Dim<2>::Tensor Tensor2d;
typedef Spheral::Dim<2>::SymTensor SymTensor2d;
typedef Spheral::Dim<2>::ThirdRankTensor ThirdRankTensor2d;
typedef Spheral::Dim<2>::FourthRankTensor FourthRankTensor2d;
typedef Spheral::Dim<2>::FifthRankTensor FifthRankTensor2d;
typedef Spheral::Dim<2>::FacetedVolume FacetedVolume2d;
typedef Spheral::GeomPlane<Spheral::Dim<3>> Plane3d;
typedef Spheral::Dim<3>::Vector Vector3d;
typedef Spheral::Dim<3>::Tensor Tensor3d;
typedef Spheral::Dim<3>::SymTensor SymTensor3d;
typedef Spheral::Dim<3>::ThirdRankTensor ThirdRankTensor3d;
typedef Spheral::Dim<3>::FourthRankTensor FourthRankTensor3d;
typedef Spheral::Dim<3>::FifthRankTensor FifthRankTensor3d;
typedef Spheral::Dim<3>::FacetedVolume FacetedVolume3d;

namespace Spheral {
namespace PYB11utils {

// Helper method for converting to/from STL vectors and Python lists
template<typename T>
inline
py::list
to_list(const std::vector<T>& stuff) {
  py::list result;
  for (const auto& x: stuff) result.append(x);
  return result;
}

template<typename T>
inline
std::vector<T>
from_list(const py::list& stuff) {
  std::vector<T> result;
  for (auto obj: stuff) result.push_back(obj.cast<T>());
  return result;
}

}
}

"""

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
               "std::vector<std::vector<std::string>>",

               "std::pair<double, double>",
               "std::pair<double, std::string>",
               "std::pair<unsigned, unsigned>",
               "std::pair<uint64_t, uint64_t>",
               "std::pair<std::string, std::string>",

               "std::map<std::string, double>",
               "std::map<int, std::string>",

               "std::vector<Dim<1>::Vector>",
               "std::vector<Dim<1>::Tensor>",
               "std::vector<Dim<1>::SymTensor>",
               "std::vector<Dim<1>::ThirdRankTensor>",
               "std::vector<Dim<1>::FourthRankTensor>",
               "std::vector<Dim<1>::FifthRankTensor>",
               "std::vector<Dim<1>::FacetedVolume>",

               "std::vector<Dim<2>::Vector>",
               "std::vector<Dim<2>::Tensor>",
               "std::vector<Dim<2>::SymTensor>",
               "std::vector<Dim<2>::ThirdRankTensor>",
               "std::vector<Dim<2>::FourthRankTensor>",
               "std::vector<Dim<2>::FifthRankTensor>",
               "std::vector<Dim<2>::FacetedVolume>",

               "std::vector<Dim<3>::Vector>",
               "std::vector<Dim<3>::Tensor>",
               "std::vector<Dim<3>::SymTensor>",
               "std::vector<Dim<3>::ThirdRankTensor>",
               "std::vector<Dim<3>::FourthRankTensor>",
               "std::vector<Dim<3>::FifthRankTensor>",
               "std::vector<Dim<3>::FacetedVolume>",

               "std::vector<GeomFacet2d>",
               "std::vector<GeomFacet3d>",

               "std::vector<Plane1d>",
               "std::vector<Plane2d>",
               "std::vector<Plane3d>",

               "std::vector<PolyClipper::Vertex2d<>>",
               "std::vector<PolyClipper::Vertex3d<>>",

               "std::vector<RKCoefficients<Dim<1>>>",
               "std::vector<RKCoefficients<Dim<2>>>",
               "std::vector<RKCoefficients<Dim<3>>>"]
