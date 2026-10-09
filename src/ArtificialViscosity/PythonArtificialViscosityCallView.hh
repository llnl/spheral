//---------------------------------Spheral++----------------------------------//
// PythonArtificialViscosityCallView -- A small view class type for the
// PythonArtificialViscosity prototyper class.
//
// Created by LDO and some agents, October 2026
//----------------------------------------------------------------------------//
#ifndef __Spheral_PythonArtificialViscosityCallView__
#define __Spheral_PythonArtificialViscosityCallView__

#include "Field/FieldList.hh"

#include <tuple>

namespace Spheral {

template<typename Dimension, typename QPiType> class PythonArtificialViscosity;

//------------------------------------------------------------------------------
// CPU-only adapter that invokes the Python-overridable viscosity method.
// It intentionally owns no viscosity parameters: those belong to its parent.
//------------------------------------------------------------------------------
template<typename Dimension, typename QPiType>
class PythonArtificialViscosityCallView {
public:
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ReturnType = QPiType;

  PythonArtificialViscosityCallView(
    const PythonArtificialViscosity<Dimension, QPiType>* parent);

  void QPiij(QPiType& QPiij, QPiType& QPiji,
             Scalar& Qij, Scalar& Qji,
             const size_t nodeListi, const size_t i,
             const size_t nodeListj, const size_t j,
             const Vector& xi,
             const SymTensor& Hi,
             const Vector& etai,
             const Vector& vi,
             const Scalar rhoi,
             const Scalar csi,
             const Vector& xj,
             const SymTensor& Hj,
             const Vector& etaj,
             const Vector& vj,
             const Scalar rhoj,
             const Scalar csj,
             const FieldListView<Dimension, Scalar>& fCl,
             const FieldListView<Dimension, Scalar>& fCq,
             const FieldListView<Dimension, Tensor>& DvDx) const;

private:
  const PythonArtificialViscosity<Dimension, QPiType>* mParent;
};

}

#endif
