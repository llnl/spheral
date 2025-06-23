//---------------------------------Spheral++----------------------------------//
// DEM -- damped linear spring contact model based on pkdgrav immplementation
//----------------------------------------------------------------------------
// Schwartz, S.R. and Richards, D.C. "An implementation of the soft-sphere 
// discrete element method in a high-performance parallel gravity tree-code,"
// Granular Matter, (2012) 14:363–380, 10.1007/s10035-012-0346-z.
//
// Zhang et. al. "Rotational Failure of Rubble-pile Bodies: Influences of 
// Shear and Cohesive Strengths," The Astrophysical Journal, (2018) 857:15, 20
// 10.3847/1538-4357/aab5b2.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#ifndef __Spheral_LinearSpringDEM_hh__
#define __Spheral_LinearSpringDEM_hh__

#include "Geometry/Dimension.hh"
#include "DEM/DEMDimension.hh"
#include "DEMBase.hh"
#include <string>

namespace Spheral {

template<typename Dimension> class DataBase;
template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
class FileIO;

template<typename Dimension>
class LinearSpringDEM : public DEMBase<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;

  typedef typename DEMBase<Dimension>::TimeStepType TimeStepType;
  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  
  LinearSpringDEM(const DataBase<Dimension>& dataBase,
                  const Scalar normalSpringConstant,
                  const Scalar normalRestitutionCoefficient,
                  const Scalar tangentialSpringConstant,
                  const Scalar tangentialRestitutionCoefficient,
                  const Scalar dynamicFrictionCoefficient,
                  const Scalar staticFrictionCoefficient,
                  const Scalar rollingFrictionCoefficient,
                  const Scalar torsionalFrictionCoefficient,
                  const Scalar cohesiveTensileStrength,
                  const Scalar shapeFactor,
                  const Scalar stepsPerCollision,
                  const bool enableFastTimeStepping,
                  const Vector& xmin,
                  const Vector& xmax);

  ~LinearSpringDEM();
 
  // physics package methods
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar time) const override;

  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                       State<Dimension>& state,
                                       StateDerivatives<Dimension>& derivs) override;

  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                         StateDerivatives<Dimension>& derivs) const override;

  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;

  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // sub-methods for dt
  TimeStepType variableTimeStep(const DataBase<Dimension>& dataBase,
                                const State<Dimension>& state,
                                const StateDerivatives<Dimension>& derivs,
                                const Scalar time) const;

  TimeStepType fixedTimeStep() const;

  void recomputeContactDuration();

  // generalized spring damper functions (inlined)
  void slidingSpringDamper(const Scalar  k,
                           const Scalar  C,
                           const Scalar  mus,
                           const Scalar  mud,
                           const Vector& x,
                           const Vector& DxDt,
                           const Scalar  fnMag,
                           const Scalar  invK,
                           const Vector& rhatij,
                           const bool    allowSiding,
                                 Vector& xNew,
                                 Vector& force) const;
                                       
  void slidingSpringDamper(const Scalar  k,
                           const Scalar  C,
                           const Scalar  mus,
                           const Scalar  mud,
                           const Scalar  x,
                           const Scalar  DxDt,
                           const Scalar  fnMag,
                           const Scalar  invK,
                           const bool    allowSiding,
                                 Scalar& xNew,
                                 Scalar& force) const;

  // set/gets
  bool enableFastTimeStepping() const;
  void   enableFastTimeStepping(bool x);

  Scalar normalSpringConstant() const;
  void   normalSpringConstant(Scalar x);

  Scalar normalRestitutionCoefficient() const;
  void   normalRestitutionCoefficient(Scalar x);

  Scalar tangentialSpringConstant() const;
  void   tangentialSpringConstant(Scalar x);

  Scalar tangentialRestitutionCoefficient() const;
  void   tangentialRestitutionCoefficient(Scalar x);

  Scalar dynamicFrictionCoefficient() const;
  void   dynamicFrictionCoefficient(Scalar x);

  Scalar staticFrictionCoefficient() const;
  void   staticFrictionCoefficient(Scalar x);

  Scalar rollingFrictionCoefficient() const;
  void   rollingFrictionCoefficient(Scalar x);

  Scalar torsionalFrictionCoefficient() const;
  void   torsionalFrictionCoefficient(Scalar x);

  Scalar cohesiveTensileStrength() const;
  void   cohesiveTensileStrength(Scalar x);

  Scalar shapeFactor() const;
  void   shapeFactor(Scalar x);

  Scalar normalBeta() const;
  void   normalBeta(Scalar x);

  Scalar tangentialBeta() const;
  void   tangentialBeta(Scalar x);

  Scalar collisionDuration() const;
  void   collisionDuration(Scalar x);

  // set moment of inertia on start up
  void setMomentOfInertia();

  // specialization for moment of Inertia different dimension
  Scalar momentOfInertia(const Scalar massi,
                         const Scalar particleRadiusi) const;

  // get methods for class FieldLists
  const FieldList<Dimension,Scalar>& momentOfInertia() const;
  const FieldList<Dimension,Scalar>& maximumOverlap() const;
  const FieldList<Dimension,Scalar>& newMaximumOverlap() const;

  
  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "LinearSpringDEM" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const override;
  virtual void restoreState(const FileIO& file, const std::string& pathName) override;
  //****************************************************************************
private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mEnableFastTimeStepping;
  Scalar mNormalSpringConstant;
  Scalar mNormalRestitutionCoefficient;
  Scalar mTangentialSpringConstant;
  Scalar mTangentialRestitutionCoefficient;
  Scalar mDynamicFrictionCoefficient;       // coefficient of friciton - dynamic
  Scalar mStaticFrictionCoefficient;        // coefficient of friction - static
  Scalar mRollingFrictionCoefficient;       // coefficient of friction - rolling
  Scalar mTorsionalFrictionCoefficient;     // coefficient of friction - torsional 
  Scalar mCohesiveTensileStrength;
  Scalar mShapeFactor;                      // varies between 0 and 1 to account to non spherical shapes & influences rolling/torsion spring parameters

  Scalar mNormalBeta;
  Scalar mTangentialBeta;
  Scalar mCollisionDuration;

  // field Lists
  FieldList<Dimension,Scalar> mMomentOfInertia;
  FieldList<Dimension,Scalar> mMaximumOverlap;
  FieldList<Dimension,Scalar> mNewMaximumOverlap;
//  FieldList<Dimension,Scalar> mOptimalSpringConstant;

  // No default constructor, copying, or assignment.
  LinearSpringDEM();
  LinearSpringDEM(const LinearSpringDEM&);
  LinearSpringDEM& operator=(const LinearSpringDEM&);
};

}

#include "LinearSpringDEMInline.hh"

#endif
