//---------------------------------Spheral++----------------------------------//
// DEMBase -- The DEM base class package for Spheral++. THis class implements
//            much of the machinery needs for DEM physics packages which 
//            should all be derived from this class or one of its daughters
//
// J.M. Pearl 2022
// ---------------------------------------------------------------------------//
// numParticleParticleContacts is currently kind of dangerous since it comes 
//     from the pair list not the contact list.
//----------------------------------------------------------------------------//

#ifndef __Spheral_DEMBase_hh__
#define __Spheral_DEMBase_hh__

#include <string>
#include "Physics/Physics.hh"
#include "DEM/DEMDimension.hh"
#include "DEM/ContactStorageLocation.hh"

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class SolidBoundaryBase;
class FileIO;
class RedistributionNotificationHandle;
//struct ContactIndex;

template<typename Dimension>
class DEMBase: public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename DEMDimension<Dimension>::AngularVector RotationType;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;
  typedef typename std::vector<SolidBoundaryBase<Dimension>*>::iterator SolidBoundaryIterator;
  typedef typename std::vector<SolidBoundaryBase<Dimension>*>::const_iterator ConstSolidBoundaryIterator;
  typedef std::shared_ptr<RedistributionNotificationHandle> RedistributionRegistrationType;
  
  // Constructors.
  DEMBase(const DataBase<Dimension>& dataBase,
          const Scalar stepsPerCollision,
          const Vector& xmin,
          const Vector& xmax);

  // Destructor.
  virtual ~DEMBase();

  // Tasks we do once on problem startup.
  virtual
  void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // Register the state Hydro expects to use and evolve.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  // This method is called once at the beginning of a timestep, after all state registration.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  // Initialize the Hydro before we start a derivative evaluation.
  virtual
  bool initialize(const Scalar time,
                  const Scalar dt,
                  const DataBase<Dimension>& dataBase,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) override;
                       

  // Finalize the derivatives.
  virtual
  void finalizeDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;
  
  // hook after the intergrator step
  // virtual
  // void finalize(const Scalar time,
  //               const Scalar dt,
  //               const DataBase<Dimension>& dataBase,
  //               const State<Dimension>& state,
  //                     StateDerivatives<Dimension>& derivs) const override;

  // Apply boundary conditions to the physics specific fields.
  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  // Enforce boundary conditions for the physics specific fields.
  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;
  
  // methods needed to square our pair fieldlists before redistribution
  void initializeBeforeRedistribution();
  void finalizeAfterRedistribution();

  void prepNeighborIndicesForRedistribution();
  
  template<typename Value>
  void prepPairFieldListForRedistribution(Value& pairFieldList);

  // templated methods
  template<typename Value1, typename Value2>
  void addContactsToPairFieldList(Value1& pairFieldList, const Value2& newValue) const;

  template<typename Value>
  void removeInactiveContactsFromPairFieldList(Value& pairFieldList) const;

  //#############################################################################
  // special methods for the pair fields if you add pairFieldLists 
  // make sure to follow the pattern in these methods.
  virtual
  void resizePairFieldLists();

  virtual
  void resizeDerivativePairFieldLists(StateDerivatives<Dimension>& derivs) const;

  virtual
  void resizeStatePairFieldLists(State<Dimension>& state) const;

  virtual 
  void removeInactiveContactsFromPairFieldLists();

  virtual 
  void removeInactiveContactsFromStatePairFieldLists(State<Dimension>& state) const;

  virtual 
  void removeInactiveContactsFromDerivativePairFieldLists(StateDerivatives<Dimension>& derivs) const;
  //#############################################################################

  void initializeOverlap(const DataBase<Dimension>& dataBase, const int startingCompositeParticleIndex);
  
  void updateContactMap(const DataBase<Dimension>& dataBase);
  
  void identifyInactiveContacts(const DataBase<Dimension>& dataBase);

  void updatePairwiseFieldLists(const bool purgeInactiveContacts = false);

  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  Scalar stepsPerCollision() const;
  void   stepsPerCollision(Scalar x);

  int  cycle() const;
  void cycle(int x);

  int  contactRemovalFrequency() const;
  void contactRemovalFrequency(int x);

  // access the dataBase 
  const DataBase<Dimension>& dataBase() const;

  // access for node fieldLists
  const FieldList<Dimension, int>&    timeStepMask() const;
  const FieldList<Dimension, Vector>& DxDt() const;
  const FieldList<Dimension, Vector>& DvDt() const;
  const FieldList<Dimension, RotationType>& omega() const;
  const FieldList<Dimension, RotationType>& DomegaDt() const;

  // access for pair fieldLists
  const FieldList<Dimension, std::vector<int>>&    isActiveContact() const;
  const FieldList<Dimension, std::vector<size_t>>& neighborIndices() const;
  const FieldList<Dimension, std::vector<Vector>>& shearDisplacement() const;
  const FieldList<Dimension, std::vector<Vector>>& rollingDisplacement() const;
  const FieldList<Dimension, std::vector<Scalar>>& torsionalDisplacement() const;
  const FieldList<Dimension, std::vector<Vector>>& DDtShearDisplacement() const;
  const FieldList<Dimension, std::vector<Vector>>& newShearDisplacement() const;
  const FieldList<Dimension, std::vector<Vector>>& DDtRollingDisplacement() const;
  const FieldList<Dimension, std::vector<Vector>>& newRollingDisplacement() const;
  const FieldList<Dimension, std::vector<Scalar>>& DDtTorsionalDisplacement() const;
  const FieldList<Dimension, std::vector<Scalar>>& newTorsionalDisplacement() const;
  const FieldList<Dimension, std::vector<Scalar>>& equilibriumOverlap() const;
  
  const std::vector<ContactIndex>& contactStorageIndices() const;


  // inlined and specialized for different dimensions
  RotationType torsionMoment(const Vector rhatij,
                             const RotationType omegai,
                             const RotationType omegaj) const;

  RotationType rollingMoment(const Vector rhatij,
                             const Vector vroti,
                             const Vector vrotj) const;

  // Solid Bounderies 
  int newSolidBoundaryIndex() const;
  void appendSolidBoundary(SolidBoundaryBase<Dimension>& boundary);
  void clearSolidBoundaries();
  void removeSolidBoundary(const SolidBoundaryBase<Dimension>& boundary);
  bool haveSolidBoundary(const SolidBoundaryBase<Dimension>& boundary) const;
  unsigned int numSolidBoundaries() const;
  const std::vector<SolidBoundaryBase<Dimension>*>& solidBoundaryConditions() const;

  // counts
  unsigned int numParticleParticleContacts() const;
  unsigned int numParticleBoundaryContacts() const;
  unsigned int numContacts() const;

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "DEMBase" ; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

protected:
  //---------------------------  Protected Interface ---------------------------//

  const DataBase<Dimension>& mDataBase;

  int mNewSolidBoundaryIndex;
  std::vector<SolidBoundaryBase<Dimension>*> mSolidBoundaries;

  int mCycle;                    // current cycle
  int mContactRemovalFrequency;  // how often do we clear out old contacts

  // number of steps per collision time-scale
  Scalar mStepsPerCollision;        

  // Optional bounding box for generating the mesh.
  Vector mxmin, mxmax;

  // fields attached to the nodes
  FieldList<Dimension, int>          mTimeStepMask;   // filter out particle from timestep calc (not active right now for DEM)
  FieldList<Dimension, Vector>       mDxDt;           // linear position derivative
  FieldList<Dimension, Vector>       mDvDt;           // linear acceleration
  FieldList<Dimension, RotationType> mOmega;          // angular velocity
  FieldList<Dimension, RotationType> mDomegaDt;       // angular acceleration

  // state fields attached to the pair interactions
  FieldList<Dimension,std::vector<size_t>> mNeighborIndices;           // tracks unique indices of contacts-we upate these (note treated specially compared to other state pair field lists)
  FieldList<Dimension,std::vector<Scalar>> mEquilibriumOverlap;        // nonzero values for composite particles
  FieldList<Dimension,std::vector<Vector>> mShearDisplacement;         // displacement for sliding spring
  FieldList<Dimension,std::vector<Vector>> mRollingDisplacement;       // displacement for rolling spring
  FieldList<Dimension,std::vector<Scalar>> mTorsionalDisplacement;     // displacement for torsional spring
  FieldList<Dimension,std::vector<int>>    mIsActiveContact;           // tracks if a interaction is still active

   // deriv fields attached to the pair interactions
  FieldList<Dimension,std::vector<Vector>> mDDtShearDisplacement;      // derivative to evolve frictional spring displacement
  FieldList<Dimension,std::vector<Vector>> mNewShearDisplacement;      // handles rotation of frictional spring and reset on slip
  FieldList<Dimension,std::vector<Vector>> mDDtRollingDisplacement;    // derivative to evolve frictional spring displacement
  FieldList<Dimension,std::vector<Vector>> mNewRollingDisplacement;    // handles rotation of frictional spring and reset on slip
  FieldList<Dimension,std::vector<Scalar>> mDDtTorsionalDisplacement;  // derivative to evolve frictional spring displacement
  FieldList<Dimension,std::vector<Scalar>> mNewTorsionalDisplacement;  // handles rotation of frictional spring and reset on slip

  // map to storage location from connectivityMap to pairwise fieldlists
  std::vector<ContactIndex> mContactStorageIndices;

  // The restart registration.
  RestartRegistrationType mRestart;
  RedistributionRegistrationType mRedistribute;

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor, copying, or assignment.
  DEMBase();
  DEMBase(const DEMBase&);
  DEMBase& operator=(const DEMBase&);
};

}

#include "DEMBaseInline.hh"

#endif
