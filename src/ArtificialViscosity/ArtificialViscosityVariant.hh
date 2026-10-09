//---------------------------------Spheral++----------------------------------//
// ArtificialViscosityVariant -- The closed set of artificial-viscosity views.
//
// This variant is visited on the host before launching a device kernel.  GPU
// builds contain only device-callable value views.  CPU-only builds additionally
// contain the Python callback adapters.
//----------------------------------------------------------------------------//
#ifndef __Spheral_ArtificialViscosityVariant__
#define __Spheral_ArtificialViscosityVariant__

#include "MonaghanGingoldViscosityView.hh"
#include "LimitedMonaghanGingoldViscosityView.hh"
#include "FiniteVolumeViscosityView.hh"
#include "TensorMonaghanGingoldViscosityView.hh"

#include <variant>

#if !defined(SPHERAL_ENABLE_HIP) && !defined(SPHERAL_ENABLE_CUDA)
#include "PythonArtificialViscosity.hh"
#endif

namespace Spheral {

template<typename Dimension>
using DeviceArtificialViscosityVariant = std::variant<
  std::monostate,
  MonaghanGingoldViscosityView<Dimension>,
  LimitedMonaghanGingoldViscosityView<Dimension>,
  FiniteVolumeViscosityView<Dimension>,
  TensorMonaghanGingoldViscosityView<Dimension>>;

#if !defined(SPHERAL_ENABLE_HIP) && !defined(SPHERAL_ENABLE_CUDA)
template<typename Dimension>
using ArtificialViscosityVariant = std::variant<
  std::monostate,
  MonaghanGingoldViscosityView<Dimension>,
  LimitedMonaghanGingoldViscosityView<Dimension>,
  FiniteVolumeViscosityView<Dimension>,
  TensorMonaghanGingoldViscosityView<Dimension>,
  PythonArtificialViscosityCallView<Dimension, typename Dimension::Scalar>,
  PythonArtificialViscosityCallView<Dimension, typename Dimension::Tensor>>;
#else
template<typename Dimension>
using ArtificialViscosityVariant = DeviceArtificialViscosityVariant<Dimension>;
#endif

}

#endif
