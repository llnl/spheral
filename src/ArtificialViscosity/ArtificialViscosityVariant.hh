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
#include "Utilities/DBC.hh"

#include <type_traits>
#include <utility>
#include <variant>

#if !defined(SPHERAL_ENABLE_HIP) && !defined(SPHERAL_ENABLE_CUDA)
#include "PythonArtificialViscosityCallView.hh"
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

template<typename Variant, typename Visitor>
decltype(auto)
visitArtificialViscosity(Variant&& Qvariant,
                         Visitor&& visitor) {
  return std::visit(
    [&](const auto& Qview) -> decltype(auto) {
      using ViewType = std::decay_t<decltype(Qview)>;

      if constexpr (std::is_same_v<ViewType, std::monostate>) {
        VERIFY2(false,
                "A TensorSVPHViscosity cannot be used as a pairwise "
                "artificial viscosity.");
      } else {
        return std::forward<Visitor>(visitor)(Qview);
      }
    },
    std::forward<Variant>(Qvariant));
}

}

#endif
