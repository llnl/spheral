//---------------------------------Spheral++----------------------------------//
// A simple form of the reducing artificial viscosity from Morris & Monaghan.
// Computes a correction scaling for the typical Monaghan & Gingold visc.
// References:
//   Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
//   Monaghan, J. J. 1992, ARA&A, 30, 543
//   Morris, J. P., & Monaghan, J. J. 1997, J. Comput. Phys., 136, 41
//
// Created by CDR, Aug 8th, 2014
//----------------------------------------------------------------------------//
#ifndef __Spheral_MorrisMonaghanReducingViscosity__
#define __Spheral_MorrisMonaghanReducingViscosity__

#include "Physics/Physics.hh"

namespace Spheral {
    
template<typename Dimension>
class MorrisMonaghanReducingViscosity: public Physics<Dimension>{
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;

  // Constructors & Destructors
  MorrisMonaghanReducingViscosity(const Scalar nhQ,
                                  const Scalar nhL,
                                  const Scalar aMin,
                                  const Scalar aMax,
                                  const Scalar negligibleSoundSpeed = 1.0e-10);
  virtual ~MorrisMonaghanReducingViscosity() = default;
    
  // No default constructor, copying, or assignment
  MorrisMonaghanReducingViscosity() = delete;
  MorrisMonaghanReducingViscosity(const MorrisMonaghanReducingViscosity&) = delete;
  MorrisMonaghanReducingViscosity& operator=(const MorrisMonaghanReducingViscosity&) const = delete;
    
  //............................................................................
  // Physics methods.
  // Increment the derivatives.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;
        
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
    
  // Apply boundary conditions to the physics specific fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;

  // Do any required one-time initializations on problem start up.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  //............................................................................
  // Restart methods.
  virtual std::string label() const override { return "MorrisMonaghanReducingViscosity"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
        
  // Access the FieldList of Reducing Viscosity multiplicative correction.
  const FieldList<Dimension, Scalar>& ClMultiplier()  const { return mClMultiplier; }
  const FieldList<Dimension, Scalar>& CqMultiplier()  const { return mCqMultiplier; }
  const FieldList<Dimension, Scalar>& DrvAlphaDtQ()   const { return mDrvAlphaDtQ; }
  const FieldList<Dimension, Scalar>& DrvAlphaDtL()   const { return mDrvAlphaDtL; }

  Scalar nhQ()                                        const { return mnhQ; }
  Scalar nhL()                                        const { return mnhL; }
  Scalar aMin()                                       const { return maMin; }
  Scalar aMax()                                       const { return maMax; }
  Scalar negligibleSoundSpeed()                       const { return mNegCs; }
    
  void nhQ(const Scalar x)                                  { mnhQ = x; }
  void nhL(const Scalar x)                                  { mnhL = x; }
  void aMin(const Scalar x)                                 { maMin = x; }
  void aMax(const Scalar x)                                 { maMax = x; }
  void negligibleSoundSpeed(const Scalar x)                 { mNegCs = x; }
    
private:
  //--------------------------- Private Interface ---------------------------//
        
  Scalar mnhQ, mnhL, maMin, maMax, mNegCs;

  FieldList<Dimension, Scalar>    mClMultiplier;
  FieldList<Dimension, Scalar>    mCqMultiplier;
  FieldList<Dimension, Scalar>    mDrvAlphaDtQ;
  FieldList<Dimension, Scalar>    mDrvAlphaDtL;

  // The restart registration.
  RestartRegistrationType mRestart;
};
    
}

#endif
