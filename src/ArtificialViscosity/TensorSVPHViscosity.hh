//---------------------------------Spheral++----------------------------------//
// A version of our tensor viscosity specialized for the SVPHFacetedHydro
// algorithm.
//
// Created by J. Michael Owen, Sat Aug 31 13:31:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorSVPHViscosity__
#define __Spheral_TensorSVPHViscosity__

#include "ArtificialViscosity.hh"

namespace Spheral {

template<typename Dimension>
class TensorSVPHViscosity: public ArtificialViscosity<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  // Constructors, destructor
  TensorSVPHViscosity(const Scalar Clinear,
                      const Scalar Cquadratic,
                      const TableKernel<Dimension>& WT,
                      const Scalar fslice);
  virtual ~TensorSVPHViscosity() = default;

  // Initialize the artificial viscosity for all FluidNodeLists in the given
  // DataBase.
  virtual bool initialize(const Scalar t,
                          const Scalar dt,
                          const DataBase<Dimension>& dataBase,
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivs) override;

  // We are abusing the normal ArtificialViscosity interface, and this normally
  // required method is a no-op for this specialization.
  // virtual void QPiij(Tensor& QPiij, Tensor& QPiji,      // result for QPi (Q/rho^2)
  //                    Scalar& Qij, Scalar& Qji,          // result for viscous pressure
  //                    const size_t nodeListi, const size_t i, 
  //                    const size_t nodeListj, const size_t j,
  //                    const Vector& xi,
  //                    const SymTensor& Hi,
  //                    const Vector& etai,
  //                    const Vector& vi,
  //                    const Scalar rhoi,
  //                    const Scalar csi,
  //                    const Vector& xj,
  //                    const SymTensor& Hj,
  //                    const Vector& etaj,
  //                    const Vector& vj,
  //                    const Scalar rhoj,
  //                    const Scalar csj,
  //                    const FieldList<Dimension, Scalar>& fCl,
  //                    const FieldList<Dimension, Scalar>& fCq,
  //                    const FieldList<Dimension, Tensor>& DvDx) const override { VERIFY2(false, "TensorSVPHViscosity ERROR: cannot call QPiij"); }

  // Access our internal state.
  Scalar fslice()                               const          { return mfslice; }
  void fslice(const Scalar x)                                  { mfslice = x; }

  const std::vector<Tensor>& DvDx()             const          { return mDvDx; }
  const std::vector<Scalar>& shearCorrection()  const          { return mShearCorrection; }
  const std::vector<Tensor>& Qface()            const          { return mQface; }

  virtual std::type_index QPiTypeIndex() const override        { return typeid(Tensor); }

  // Restart methods.
  virtual std::string label()                   const override { return "TensorSVPHViscosity"; }

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mfslice;
  std::vector<Tensor> mDvDx;
  std::vector<Scalar> mShearCorrection;
  std::vector<Tensor> mQface;
};

}

#endif
