//---------------------------------Spheral++----------------------------------//
// NBodyGravity -- A simple implementation of the N-body gravity algorithm.
//
//! \author $Author: mikeowen $
//! \version $Revision: 3174 $
//! \date $Date: 2009-07-22 14:23:24 -0700 (Wed, 22 Jul 2009) $
//
//----------------------------------------------------------------------------//
#ifndef NBODYGRAVITY_HH
#define NBODYGRAVITY_HH

#include "Physics/GenericBodyForce.hh"
#include "Field/FieldList.hh"

#include <vector>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

template <typename Dimension>
class NBodyGravity: public GenericBodyForce<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
  using ConnectivityRequirements = typename Physics<Dimension>::ConnectivityRequirements;

  //! Constructor.
  //! \param plummerSofteningLength -- The Plummer Softening Length for the model.
  //! \param maxDeltaVelocity -- Maximum factor by which the velocity can be changed by an acceleration per timestep.
  NBodyGravity(const double plummerSofteningLength,
               const double maxDeltaVelocity,
               const double G,
               const bool compatibleVelocityUpdate);

  //! Destructor.
  virtual ~NBodyGravity();

  //! We augment the generic body force state.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  //! This is the derivative method that all BodyForce classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const override;

  //! Vote on the timestep.  This uses a velocity-limiting rule.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // An optional hook to initialize once when the problem is starting up.
  // Typically this is used to size arrays once all the materials and NodeLists have
  // been created.  It is assumed after this method has been called it is safe to
  // call Physics::registerState for instance to create full populated State objects.
  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;
  //! Beginning of timestep work.
  virtual void preStepInitialize(const DataBase<Dimension>& dataBase, 
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  //! End of timestep finalizations.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs) override;

  //! Required label for Physics interface.
  virtual std::string label() const override { return "NBodyGravity"; }

  //! This package opts out of building connectivity.
  //! Returns {connectivity, ghostConnectivity, overlapConnectivity, intersectionConnectivity}
  virtual ConnectivityRequirements requireConnectivity() const override { return {false, false, false, false}; }

  //! Return the total energy contribution due to the gravitational potential.
  virtual Scalar extraEnergy() const override;

  //! Return the gravitational potential created by the particle distribution.
  const FieldList<Dimension, Scalar>& potential() const;

  //! Test if the package is valid, i.e., ready to use.
  virtual bool valid() const;

  //! The gravitational constant we're using.
  double G() const;

  //! The current softening length.
  double softeningLength() const;
  void softeningLength(double x);

  //! Flag for using the compatible velocity update.
  bool compatibleVelocityUpdate() const;
  void compatibleVelocityUpdate(bool x);

private:
  
  //! The gravitational potential of the particles.
  mutable FieldList<Dimension, Scalar> mPotential;
  mutable FieldList<Dimension, Scalar> mPotential0;
  mutable FieldList<Dimension, Scalar> mVel02;

  //! The total potential energy of the particles.  Also mutable.
  mutable Scalar mExtraEnergy;

  //! The maximum acceleration magnitude during the last time step.
  mutable Scalar mOldMaxAcceleration;
  
  //! The maximum velocity magnitude during the last time step.
  mutable Scalar mOldMaxVelocity;

  //! The maximum allowed change in velocity, by factors of velocity.
  Scalar mMaxDeltaVelocityFactor;

  //! The Plummer softening length.
  Scalar mSofteningLength;
  
  //! The gravitational constant.
  Scalar mG;
  
  //! Flag for compatible velocity.
  bool mCompatibleVelocityUpdate;

  // Default constructor -- disabled.
  NBodyGravity();

  // Copy constructor -- disabled.
  NBodyGravity(const NBodyGravity&);

  // Assignment operator -- disabled.
  NBodyGravity& operator=(const NBodyGravity&);

  // Worker for accumulating pair-wise forces.
  void applyPairForces(const std::vector<Scalar>& otherMass,
                       const std::vector<Vector>& otherPosition,
                       const FieldList<Dimension, Vector>& position,
                       FieldList<Dimension, Vector>& DvDt,
                       FieldList<Dimension, Scalar>& potential) const;

  // Methods for serializing/deserializing point values.
  void serialize(const FieldList<Dimension, typename Dimension::Scalar>& mass,
                 const FieldList<Dimension, typename Dimension::Vector>& position,
                 std::vector<char>& buffer) const;
  void deserialize(const std::vector<char>& buffer,
                   std::vector<typename Dimension::Scalar>& mass,
                   std::vector<typename Dimension::Vector>& position) const;

}; // end class NBodyGravity

} // end namespace Spheral

#endif
