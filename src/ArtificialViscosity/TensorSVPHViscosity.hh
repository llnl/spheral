//---------------------------------Spheral++----------------------------------//
// A version of our tensor viscosity specialized for the SVPHFacetedHydro
// algorithm. This class is not an ArtificialViscosityVariant class but provides
// an unrelated variant to fulfill the virtual function requirement.
//
// Created by J. Michael Owen, Sat Aug 31 13:31:51 PDT 2013
//----------------------------------------------------------------------------//
#ifndef __Spheral_TensorSVPHViscosity__
#define __Spheral_TensorSVPHViscosity__

#include "ArtificialViscosity.hh"
#include "ArtificialViscosityVariant.hh"

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

  // Forward ArtificialViscosity's virtual parameter interface to the view.
  virtual Scalar Cl() const override { return ArtificialViscosityView<Dimension>::Cl(); }
  virtual Scalar Cq() const override { return ArtificialViscosityView<Dimension>::Cq(); }
  virtual bool balsaraShearCorrection() const override {
    return ArtificialViscosityView<Dimension>::balsaraShearCorrection();
  }
  virtual Scalar epsilon2() const override { return ArtificialViscosityView<Dimension>::epsilon2(); }
  virtual Scalar negligibleSoundSpeed() const override {
    return ArtificialViscosityView<Dimension>::negligibleSoundSpeed();
  }

  virtual void Cl(const Scalar x) override { ArtificialViscosityView<Dimension>::Cl(x); }
  virtual void Cq(const Scalar x) override { ArtificialViscosityView<Dimension>::Cq(x); }
  virtual void balsaraShearCorrection(const bool x) override {
    ArtificialViscosityView<Dimension>::balsaraShearCorrection(x);
  }
  virtual void epsilon2(const Scalar x) override { ArtificialViscosityView<Dimension>::epsilon2(x); }
  virtual void negligibleSoundSpeed(const Scalar x) override {
    ArtificialViscosityView<Dimension>::negligibleSoundSpeed(x);
  }

  // Access our internal state.
  Scalar fslice()                               const          { return mfslice; }
  void fslice(const Scalar x)                                  { mfslice = x; }

  const std::vector<Tensor>& DvDx()             const          { return mDvDx; }
  const std::vector<Scalar>& shearCorrection()  const          { return mShearCorrection; }
  const std::vector<Tensor>& Qface()            const          { return mQface; }

  // Restart methods.
  virtual std::string label()                   const override { return "TensorSVPHViscosity"; }

  // Dummy function to fulfill the variantView override
  ArtificialViscosityVariant<Dimension> variantView() const override {
    VERIFY2(false, "Cannot call variantView with TensorSVPHViscosity");
    return std::monostate;
  }

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mfslice;
  std::vector<Tensor> mDvDx;
  std::vector<Scalar> mShearCorrection;
  std::vector<Tensor> mQface;
};

}

#endif
