//---------------------------------Spheral++----------------------------------//
// Dim -- trait class for dimensional type information.
//
// Dim also carries some useful static methods, which (for nu dimensions) are
//   rootnu(x) = x**(1/nu)
//   pownu(x) = x**nu
//   pownu1(x) = x**(nu - 1)
//   pownu12(x) = x**((nu - 1)/2)
//
// Created by JMO, Fri Apr 16 10:13:52 PDT 1999
//----------------------------------------------------------------------------//

#ifndef __Spheral_Dim_hh__
#define __Spheral_Dim_hh__

#include "GeomVector.hh"
#include "GeomTensor.hh"
#include "GeomSymmetricTensor.hh"
#include "Geom3Vector.hh"
#include "GeomThirdRankTensor.hh"
#include "GeomFourthRankTensor.hh"
#include "GeomFifthRankTensor.hh"
#include "Box1d.hh"
#include "GeomPolygon.hh"
#include "GeomPolyhedron.hh"
#include "Utilities/FastMath.hh"

#include <cmath>

namespace Spheral {

template<int ndim>
class Dim {
};

template<>
class Dim<1> {
public:
  typedef double Scalar;
  typedef GeomVector<1> Vector;
  typedef Geom3Vector Vector3d;
  typedef GeomTensor<1> Tensor;
  typedef GeomSymmetricTensor<1> SymTensor;
  typedef GeomThirdRankTensor<1> ThirdRankTensor;
  typedef GeomFourthRankTensor<1> FourthRankTensor;
  typedef GeomFifthRankTensor<1> FifthRankTensor;
  typedef GeomVector<1> WMVector;
  typedef Box1d Box;
  typedef Box1d ConvexHull;
  typedef Box1d FacetedVolume;
  typedef GeomFacet1d Facet;
  static const int nDim = 1;

  static SPHERAL_HOST_DEVICE double rootnu(const double& x) { return x; }
  static double pownu(const double& x) { return x; }
  static double pownu1(const double& /*x*/) { return 1.0; }
  static double pownu12(const double& /*x*/) { return 1.0; }
};

template<>
class Dim<2> {
public:
  typedef double Scalar;
  typedef GeomVector<2> Vector;
  typedef Geom3Vector Vector3d;
  typedef GeomTensor<2> Tensor;
  typedef GeomSymmetricTensor<2> SymTensor;
  typedef GeomThirdRankTensor<2> ThirdRankTensor;
  typedef GeomFourthRankTensor<2> FourthRankTensor;
  typedef GeomFifthRankTensor<2> FifthRankTensor;
  typedef GeomPolygon ConvexHull;
  typedef GeomPolygon FacetedVolume;
  typedef GeomFacet2d Facet;
  static const int nDim = 2;

  static SPHERAL_HOST_DEVICE double rootnu(const double& x) { return std::sqrt(x); }
  static double pownu(const double& x) { return x*x; }
  static double pownu1(const double& x) { return x; }
  static double pownu12(const double& x) { return std::sqrt(x); }
};

template<>
class Dim<3> {
public:
  typedef double Scalar;
  typedef GeomVector<3> Vector;
  typedef Geom3Vector Vector3d;
  typedef GeomTensor<3> Tensor;
  typedef GeomSymmetricTensor<3> SymTensor;
  typedef GeomThirdRankTensor<3> ThirdRankTensor;
  typedef GeomFourthRankTensor<3> FourthRankTensor;
  typedef GeomFifthRankTensor<3> FifthRankTensor;
  typedef GeomPolyhedron ConvexHull;
  typedef GeomPolyhedron FacetedVolume;
  typedef GeomFacet3d Facet;
  static const int nDim = 3;

  static SPHERAL_HOST_DEVICE double rootnu(const double& x) { return FastMath::CubeRootHalley2(x); } // { return pow(x, 1.0/3.0); }
  static double pownu(const double& x) { return x*x*x; }
  static double pownu1(const double& x) { return x*x; }
  static double pownu12(const double& x) { return x; }
};

//--------------------------------------------------------------------------------
template<typename DataType1, typename DataType2>
class CombineTypes {};

template<>
class CombineTypes<Dim<1>::Scalar, Dim<1>::Scalar> {
public:
  typedef Dim<1>::Scalar ProductType;
};

template<>
class CombineTypes<Dim<1>::Scalar, Dim<1>::Vector> {
public:
  typedef Dim<1>::Vector ProductType;
};

template<>
class CombineTypes<Dim<1>::Scalar, Dim<1>::Tensor> {
public:
  typedef Dim<1>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<1>::Scalar, Dim<1>::SymTensor> {
public:
  typedef Dim<1>::SymTensor ProductType;
};

template<>
class CombineTypes<Dim<1>::Vector, Dim<1>::Scalar> {
public:
  typedef Dim<1>::Vector ProductType;
};

template<>
class CombineTypes<Dim<1>::Vector, Dim<1>::Vector> {
public:
  typedef Dim<1>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<1>::Vector, Dim<1>::Tensor> {
public:
  typedef Dim<1>::Vector ProductType;
};

template<>
class CombineTypes<Dim<1>::Vector, Dim<1>::SymTensor> {
public:
  typedef Dim<1>::Vector ProductType;
};

template<>
class CombineTypes<Dim<1>::Tensor, Dim<1>::Scalar> {
public:
  typedef Dim<1>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<1>::Tensor, Dim<1>::Vector> {
public:
  typedef Dim<1>::Vector ProductType;
};

template<>
class CombineTypes<Dim<1>::Tensor, Dim<1>::Tensor> {
public:
  typedef Dim<1>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<1>::Tensor, Dim<1>::SymTensor> {
public:
  typedef Dim<1>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<1>::SymTensor, Dim<1>::Scalar> {
public:
  typedef Dim<1>::SymTensor ProductType;
};

template<>
class CombineTypes<Dim<1>::SymTensor, Dim<1>::Vector> {
public:
  typedef Dim<1>::Vector ProductType;
};

template<>
class CombineTypes<Dim<1>::SymTensor, Dim<1>::Tensor> {
public:
  typedef Dim<1>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<1>::SymTensor, Dim<1>::SymTensor> {
public:
  typedef Dim<1>::Tensor ProductType;
};

//--------------------------------------------------------------------------------
template<>
class CombineTypes<Dim<2>::Scalar, Dim<2>::Vector> {
public:
  typedef Dim<2>::Vector ProductType;
};

template<>
class CombineTypes<Dim<2>::Scalar, Dim<2>::Tensor> {
public:
  typedef Dim<2>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<2>::Scalar, Dim<2>::SymTensor> {
public:
  typedef Dim<2>::SymTensor ProductType;
};

template<>
class CombineTypes<Dim<2>::Vector, Dim<2>::Scalar> {
public:
  typedef Dim<2>::Vector ProductType;
};

template<>
class CombineTypes<Dim<2>::Vector, Dim<2>::Vector> {
public:
  typedef Dim<2>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<2>::Vector, Dim<2>::Tensor> {
public:
  typedef Dim<2>::Vector ProductType;
};

template<>
class CombineTypes<Dim<2>::Vector, Dim<2>::SymTensor> {
public:
  typedef Dim<2>::Vector ProductType;
};

template<>
class CombineTypes<Dim<2>::Tensor, Dim<2>::Scalar> {
public:
  typedef Dim<2>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<2>::Tensor, Dim<2>::Vector> {
public:
  typedef Dim<2>::Vector ProductType;
};

template<>
class CombineTypes<Dim<2>::Tensor, Dim<2>::Tensor> {
public:
  typedef Dim<2>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<2>::Tensor, Dim<2>::SymTensor> {
public:
  typedef Dim<2>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<2>::SymTensor, Dim<2>::Scalar> {
public:
  typedef Dim<2>::SymTensor ProductType;
};

template<>
class CombineTypes<Dim<2>::SymTensor, Dim<2>::Vector> {
public:
  typedef Dim<2>::Vector ProductType;
};

template<>
class CombineTypes<Dim<2>::SymTensor, Dim<2>::Tensor> {
public:
  typedef Dim<2>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<2>::SymTensor, Dim<2>::SymTensor> {
public:
  typedef Dim<2>::Tensor ProductType;
};

//--------------------------------------------------------------------------------
template<>
class CombineTypes<Dim<3>::Scalar, Dim<3>::Vector> {
public:
  typedef Dim<3>::Vector ProductType;
};

template<>
class CombineTypes<Dim<3>::Scalar, Dim<3>::Tensor> {
public:
  typedef Dim<3>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<3>::Scalar, Dim<3>::SymTensor> {
public:
  typedef Dim<3>::SymTensor ProductType;
};

template<>
class CombineTypes<Dim<3>::Vector, Dim<3>::Scalar> {
public:
  typedef Dim<3>::Vector ProductType;
};

template<>
class CombineTypes<Dim<3>::Vector, Dim<3>::Vector> {
public:
  typedef Dim<3>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<3>::Vector, Dim<3>::Tensor> {
public:
  typedef Dim<3>::Vector ProductType;
};

template<>
class CombineTypes<Dim<3>::Vector, Dim<3>::SymTensor> {
public:
  typedef Dim<3>::Vector ProductType;
};

template<>
class CombineTypes<Dim<3>::Tensor, Dim<3>::Scalar> {
public:
  typedef Dim<3>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<3>::Tensor, Dim<3>::Vector> {
public:
  typedef Dim<3>::Vector ProductType;
};

template<>
class CombineTypes<Dim<3>::Tensor, Dim<3>::Tensor> {
public:
  typedef Dim<3>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<3>::Tensor, Dim<3>::SymTensor> {
public:
  typedef Dim<3>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<3>::SymTensor, Dim<3>::Scalar> {
public:
  typedef Dim<3>::SymTensor ProductType;
};

template<>
class CombineTypes<Dim<3>::SymTensor, Dim<3>::Vector> {
public:
  typedef Dim<3>::Vector ProductType;
};

template<>
class CombineTypes<Dim<3>::SymTensor, Dim<3>::Tensor> {
public:
  typedef Dim<3>::Tensor ProductType;
};

template<>
class CombineTypes<Dim<3>::SymTensor, Dim<3>::SymTensor> {
public:
  typedef Dim<3>::Tensor ProductType;
};

}

#endif
