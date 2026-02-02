//---------------------------------Spheral++----------------------------------//
// An implementation of the Cullen Dehnen Viscosity, with Hopkins alterations
// Computes a correction for the scaling.
// References:
// Cullen, L., & Dehnen, W. 2010, MNRAS, 408, 669
// Hopkins arXiv:1409.7395
//----------------------------------------------------------------------------//
#ifndef __Spheral_CullenDehnenViscosity__
#define __Spheral_CullenDehnenViscosity__

#include "Physics/Physics.hh"

namespace Spheral {
    
// Forward declarations
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class SmoothingScaleBase;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
class FileIO;

template<typename Dimension>
class CullenDehnenViscosity: public Physics<Dimension>{
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;
        
  // Constructors & Destructors
  CullenDehnenViscosity(const TableKernel<Dimension>& WT,
                        const Scalar alphMax,
                        const Scalar alphMin,
                        const Scalar betaC,
                        const Scalar betaD,
                        const Scalar betaE,
                        const Scalar fKern,
                        const bool boolHopkins);
  virtual ~CullenDehnenViscosity() = default;
    
  // No default constructor, copying, or assignment
  CullenDehnenViscosity() = delete;
  CullenDehnenViscosity(const CullenDehnenViscosity&) = delete;
  CullenDehnenViscosity& operator=(const CullenDehnenViscosity&) const = delete;

  //............................................................................
  // Physics methods.
  // Increment the derivatives.
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivatives) const override;
  virtual
  void finalize(const Scalar time,
                const Scalar dt,
                DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivs) override;
  virtual   
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;
    
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
  //............................................................................
    
  // Restart methods.
  virtual std::string label() const override { return "CullenDehnenViscosity"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);

  // Access the stored data and parameters
  const TableKernel<Dimension>& kernel()               const { return mWT; }

  Scalar alphMax()                                     const { return malphMax; }
  Scalar alphMin()                                     const { return malphMin; }
  Scalar betaE()                                       const { return mbetaE; }
  Scalar betaD()                                       const { return mbetaD; }
  Scalar betaC()                                       const { return mbetaC; }
  Scalar fKern()                                       const { return mfKern; }
  bool boolHopkins()                                   const { return mboolHopkins; }
    
  void alphMax(Scalar x)                                     { malphMax = x; }
  void alphMin(Scalar x)                                     { malphMin = x; }
  void betaE(Scalar x)                                       { mbetaE = x; }
  void betaD(Scalar x)                                       { mbetaD = x; }
  void betaC(Scalar x)                                       { mbetaC = x; }
  void fKern(Scalar x)                                       { mfKern = x; }
  void boolHopkins(bool x)                                   { mboolHopkins = x; }

  const FieldList<Dimension, Scalar>& ClMultiplier()   const { return mClMultiplier; }
  const FieldList<Dimension, Scalar>& CqMultiplier()   const { return mCqMultiplier; }
  const FieldList<Dimension, Vector>& PrevDvDt()       const { return mPrevDvDt; }
  const FieldList<Dimension, Scalar>& PrevDivV()       const { return mPrevDivV; }
  const FieldList<Dimension, Scalar>& CullAlpha()      const { return mCullAlpha; }
  const FieldList<Dimension, Scalar>& PrevDivV2()      const { return mPrevDivV2; }
  const FieldList<Dimension, Scalar>& CullAlpha2()     const { return mCullAlpha2; }
  const FieldList<Dimension, Scalar>& DalphaDt()       const { return mDalphaDt; }
  const FieldList<Dimension, Scalar>& alphaLocal()     const { return mAlphaLocal; }
  const FieldList<Dimension, Scalar>& alpha0()         const { return mAlpha0; }
  const FieldList<Dimension, Scalar>& R()              const { return mR; }
  const FieldList<Dimension, Scalar>& vsig()           const { return mVsig; }

private:
  //--------------------------- Private Interface ---------------------------//

  const TableKernel<Dimension>&   mWT;

  FieldList<Dimension, Scalar>    mClMultiplier;
  FieldList<Dimension, Scalar>    mCqMultiplier;

  FieldList<Dimension, Vector>    mPrevDvDt;     //Will enroll as state fields
  FieldList<Dimension, Scalar>    mPrevDivV;
  FieldList<Dimension, Scalar>    mCullAlpha;
  FieldList<Dimension, Scalar>    mPrevDivV2;    //Will enroll as derivative fields so that we can store the previous time step value.
  FieldList<Dimension, Scalar>    mCullAlpha2;

  FieldList<Dimension, Scalar>    mDalphaDt;     // Time derivative of alpha
  FieldList<Dimension, Scalar>    mAlphaLocal;   // Alpha local to be filled in derivatives
  FieldList<Dimension, Scalar>    mAlpha0;       // The Hopkins form actually evolves alpha0
  
  FieldList<Dimension, Scalar>    mR;            //Will enroll as derivs fields
  FieldList<Dimension, Scalar>    mVsig;         //Will enroll as derivs fields

  Scalar malphMax, malphMin, mbetaC, mbetaD, mbetaE, mfKern;
  bool mboolHopkins;//Use Hopkins Reformulation

  // The restart registration.
  RestartRegistrationType mRestart;
};
    
}

#endif
