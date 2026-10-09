//---------------------------------Spheral++----------------------------------//
// A version of our tensor viscosity specialized for the SVPHFacetedHydro
// algorithm. This class is not an ArtificialViscosityVariant class.
// It does not provide a QPiij or a view.
//
// Created by J. Michael Owen, Sat Aug 31 13:31:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorSVPHViscosity__
#define __Spheral_TensorSVPHViscosity__

#include "ArtificialViscosity.hh"
#include "ArtificialViscosityView.hh"

namespace Spheral {

template<typename Dimension>
class TensorSVPHViscosity:
    public ArtificialViscosity<Dimension>,
    public ArtificialViscosityView<Dimension> {
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

  // Access our internal state.
  Scalar fslice()                               const          { return mfslice; }
  void fslice(const Scalar x)                                  { mfslice = x; }

  const std::vector<Tensor>& DvDx()             const          { return mDvDx; }
  const std::vector<Scalar>& shearCorrection()  const          { return mShearCorrection; }
  const std::vector<Tensor>& Qface()            const          { return mQface; }

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
