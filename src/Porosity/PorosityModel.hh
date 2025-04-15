//---------------------------------Spheral++----------------------------------//
// PorosityModel
// Base class for PorosityModels for common functionality.
//
// Created by JMO, Mon Nov 13 09:28:26 PST 2023
//----------------------------------------------------------------------------//
#ifndef __Spheral_PorosityModel__
#define __Spheral_PorosityModel__

#include "Physics/Physics.hh"
#include "NodeList/SolidNodeList.hh"
#include "DataOutput/registerWithRestart.hh"

#include <memory>

namespace Spheral {

template<typename Dimension>
class PorosityModel: 
    public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ResidualType = typename Physics<Dimension>::ResidualType;

  // Constructors, destructors.
  PorosityModel(const SolidNodeList<Dimension>& nodeList,        // The NodeList we're going apply to
                const double phi0,                               // Initial porosity
                const double cS0,                                // Reference sound speed at full density
                const double c0,                                 // Reference sound speed at initial porosity
                const double rhoS0,                              // Reference solid density
                const bool jutziStateUpdate);                    // Apply state update rules from Jutzi 2008

  PorosityModel(const SolidNodeList<Dimension>& nodeList,        // The NodeList we're going apply to
                const Field<Dimension, Scalar>& phi0,            // Initial porosity
                const double cS0,                                // Reference sound speed at full density
                const Field<Dimension, Scalar>& c0,              // Reference sound speed at initial porosity
                const double rhoS0,                              // Reference solid density
                const bool jutziStateUpdate);                    // Apply state update rules from Jutzi 2008

  virtual ~PorosityModel();

  //............................................................................
  // Physics methods.
  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // Register our state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Do any required one-time initializations on problem start up.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Return the maximum state change we care about for checking for convergence in the implicit integration methods.
  // For porosity default to just relying on the hydro to handle this, so return no vote.
  virtual ResidualType maxResidual(const DataBase<Dimension>& dataBase, 
                                   const State<Dimension>& state1,
                                   const State<Dimension>& state0,
                                   const Scalar tol) const override;
  //............................................................................

  //............................................................................
  // Methods required for restarting.
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //............................................................................

  // Access the material parameters.
  bool jutziStateUpdate() const;
  double rhoS0() const;
  double cS0() const;
  double KS0() const;
  double fdt() const;
  double maxAbsDalphaDt() const;
  const SolidNodeList<Dimension>& nodeList() const;
  const Field<Dimension, Scalar>& alpha0() const;
  const Field<Dimension, Scalar>& alpha() const;
  const Field<Dimension, Scalar>& DalphaDt() const;
  const Field<Dimension, Scalar>& solidMassDensity() const;
  const Field<Dimension, Scalar>& c0() const;

  // The optional scaling term for deviatoric stress update from Jutzi 2008.
  // Note: if jutziStateUpdate is false, these fields are not allocated.
  const Field<Dimension, Scalar>& fDS() const;
  const Field<Dimension, Scalar>& fDSnew() const;

  // Timestep multiplier
  void fdt(const double x);

  // Provide the porosity (phi) computed from the internally stored distention alpha
  Field<Dimension, Scalar> phi() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  bool mJutziStateUpdate;
  double mRhoS0, mcS0, mKS0, mfdt;
  mutable double mMaxAbsDalphaDt;
  const SolidNodeList<Dimension>& mNodeList;
  Field<Dimension, Scalar> mAlpha0, mAlpha, mDalphaDt, mSolidMassDensity, mc0, mfDS, mfDSnew;

  // The restart registration.
  RestartRegistrationType mRestart;

  // Disallow default constructor
  PorosityModel();
};

}

#include "PorosityModelInline.hh"

#endif
