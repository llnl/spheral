//---------------------------------Spheral++----------------------------------//
// PythonArtificialViscosity -- Implementation
//----------------------------------------------------------------------------//
#include "ArtificialViscosity/PythonArtificialViscosity.hh"
#include "Utilities/SpheralMessage.hh"

#include <stdexcept>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension, typename QPiType>
PythonArtificialViscosity<Dimension, QPiType>::
PythonArtificialViscosity(const Scalar Clinear,
                                const Scalar Cquadratic,
                                const TableKernel<Dimension>& kernel) :
  ArtificialViscosity<Dimension>(Clinear, Cquadratic, kernel),
  mView(nullptr) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension, typename QPiType>
PythonArtificialViscosity<Dimension, QPiType>::
~PythonArtificialViscosity() {
  if (mView) {
    mView.free();
  }
}

//------------------------------------------------------------------------------
// Get scalar view - create wrapper on demand
//------------------------------------------------------------------------------
template<typename Dimension, typename QPiType>
chai::managed_ptr<typename PythonArtificialViscosity<Dimension, QPiType>::ArtViscViewScalar>
PythonArtificialViscosity<Dimension, QPiType>::
getScalarView() {
  if constexpr (!std::is_same_v<QPiType, Scalar>) {
    return chai::managed_ptr<ArtViscViewScalar>();
  } else {
    if (!mView) {
      mView = chai::make_managed<PythonAVView>(this);
    }
    return chai::dynamic_pointer_cast<ArtViscViewScalar>(mView);
  }
}

//------------------------------------------------------------------------------
// Get tensor view - create wrapper on demand.
//------------------------------------------------------------------------------
template<typename Dimension, typename QPiType>
chai::managed_ptr<typename PythonArtificialViscosity<Dimension, QPiType>::ArtViscViewTensor>
PythonArtificialViscosity<Dimension, QPiType>::
getTensorView() {
  if constexpr (!std::is_same_v<QPiType, Tensor>) {
    return chai::managed_ptr<ArtViscViewTensor>();
  } else {
    if (!mView) {
      mView = chai::make_managed<PythonAVView>(this);
    }
    return chai::dynamic_pointer_cast<ArtViscViewTensor>(mView);
  }
}

//------------------------------------------------------------------------------
// PythonAVView Constructor
//------------------------------------------------------------------------------
template<typename Dimension, typename QPiType>
PythonArtificialViscosity<Dimension, QPiType>::PythonAVView::
PythonAVView(PythonArtificialViscosity<Dimension, QPiType>* parent) :
  ArtificialViscosityView<Dimension, QPiType>(parent->Cl(),
                                             parent->Cq(),
                                             parent->balsaraShearCorrection(),
                                             parent->epsilon2(),
                                             parent->negligibleSoundSpeed()),
  mParent(parent) {
  REQUIRE(parent != nullptr);
}

//------------------------------------------------------------------------------
// QPiij implementation - extract FieldListView data and call Python method
//------------------------------------------------------------------------------
template<typename Dimension, typename QPiType>
void
PythonArtificialViscosity<Dimension, QPiType>::PythonAVView::
QPiij(QPiType& QPiij, QPiType& QPiji,
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
      const FieldListView<Dimension, Tensor>& DvDx) const {

  // Extract multipliers from FieldListView
  // If FieldList is empty, use default multiplier of 1.0
  const auto fCli = (fCl.size() > 0) ? fCl(nodeListi, i) : 1.0;
  const auto fCqi = (fCq.size() > 0) ? fCq(nodeListi, i) : 1.0;
  const auto fClj = (fCl.size() > 0) ? fCl(nodeListj, j) : 1.0;
  const auto fCqj = (fCq.size() > 0) ? fCq(nodeListj, j) : 1.0;

  // Call simplified Python-overridable method with exception handling
  // This catches any Python exceptions and converts them to C++ errors
  try {
    std::tie(QPiij, QPiji, Qij, Qji) = mParent->computeQPiij(xi, vi, rhoi, csi,
                                                             xj, vj, rhoj, csj,
                                                             etai, etaj,
                                                             fCli, fCqi, fClj, fCqj);
  } catch (const std::exception& e) {
    // Python raised an exception - report it clearly
    VERIFY2(false, "Python computeQPiij raised exception: " << e.what());
  } catch (...) {
    // Unknown exception from Python
    VERIFY2(false, "Python computeQPiij raised unknown exception");
  }
}

}
