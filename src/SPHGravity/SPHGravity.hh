//---------------------------------Spheral++----------------------------------//
// SPHGravity -- A simple implementation of the N-body gravity algorithm.
//
//! \author $Author: jeffjohnson $
//! \version $Revision: 2239 $
//! \date $Date: 2007-05-28 23:58:39 -0700 (Mon, 28 May 2007) $
//
//----------------------------------------------------------------------------//
#ifndef SPHGravity_HH
#define SPHGravity_HH

#include "Python.h"
#include "Physics/Physics.hh"
#include "Field/FieldList.hh"
#include "Kernel/TableKernel.hh"

#include "petsc.h"
#include "petscksp.h"
#include <map>

namespace Spheral {

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;

namespace GravitySpace {

template <typename Dimension>
class SPHGravity: public Physics<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  typedef typename Physics<Dimension>::TimeStepType TimeStepType;
  typedef typename Physics<Dimension>::ConstBoundaryIterator ConstBoundaryIterator;

  //! Constructor.
  //! \param W -- An SPHKernel used to compute SPH interpolants.
  //! \param G -- Cavendish's gravitational constant expressed in the desired units.
  //! \param maxDeltaVelocity -- Maximum factor by which the velocity can be changed by an acceleration per timestep.
  //! \param safetyFactor -- The safety factor to apply to the timestep.  Must be between 0 and 1.
  SPHGravity(const TableKernel<Dimension>& W,
             Scalar G, 
             Scalar maxDeltaVelocity = 2.0,
             Scalar safetyFactor = 0.5);

  //! Destructor.
  virtual ~SPHGravity();

  //! This is the derivative method that all Physics classes must provide.
  virtual 
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                           StateDerivatives<Dimension>& derivs) const;

  //! Register the state for gravitational acceleration.
  virtual 
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state);

  // Register the derivatives/change fields for updating state.
  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs);

  // Apply boundary conditions to the fields.
  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs);

  // Enforce boundary conditions for the fields.
  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs);

  //! Vote on the timestep.  This uses a velocity-limiting rule.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const;

  //! Make sure that Gadget's internal state is initialized before cycling.
  virtual bool initialize(const Scalar& time, 
                          const Scalar& dt,
                          const DataBase<Dimension>& db, 
                          State<Dimension>& state,
                          StateDerivatives<Dimension>& derivs);

  //! Return the norm of the residual vector after the most recently converged
  //! linear solve.
  double residualNorm() const { return mResNorm; }

  //! Return the number of iterations taken in the most recent linear solve.
  int numIterations() const { return mNumIters; }

  //! Return the total energy contribution due to the gravitational potential.
  virtual Scalar extraEnergy() const;

  //! Return the gravitational potential created by the particle distribution.
  const FieldList<Dimension, Scalar>& potential() const;

  //! Test if the package is valid, i.e., ready to use.
  virtual bool valid() const;

  //! The gravitational constant we're using.
  double G() const { return mG; }

private:
  
  //! Information about the most recent linear solve.
  mutable double mResNorm;
  mutable int mNumIters;

  //! The gravitational constant.
  Scalar mG;

  //! SPH kernel.
  const TableKernel<Dimension>& mKernel;
  
  //! The maximum allowed change in velocity, by factors of velocity.
  Scalar mMaxVChangeFactor;

  //! The gravitational potential of the particles.
  mutable FieldList<Dimension, Scalar> mPotential;

  //! The total potential energy of the particles.  Also mutable.
  mutable Scalar mExtraEnergy;

  //! The minimum velocity / acceleration ratio. 
  mutable Scalar mMinViOverAi;

  //! The timestep safety factor.
  mutable Scalar mSafetyFactor;

  //! The minimum gravitational dynamic time scale.
  mutable Scalar mMinDynTimeScale;

  //! Vector factory.
  mutable PyObject* mVecFactory;

  //! Matrix factory.
  mutable PyObject* mMatFactory;

  //! Laplacian matrix. 
  mutable PyObject* mMatrix;

  //! Right-hand-side vector.
  mutable PyObject* mRHS;

  //! Map of boundary nodes that are related to periodic or distributed 
  //! boundary conditions.
  mutable std::map<std::pair<int, int>, int> mIsPeriodicOrDistributedNode;

  //! PETSc linear solver.
  mutable KSP mSolver;

  //! A global indexing scheme for all the nodes in our problem of interest.
  mutable FieldList<Dimension, int> mNodeIndices;

  //! A mapping of nodes to other nodes that overlap them.
  mutable std::map<int, int> mOverlapNodes;

  //! Compute the Laplacian matrix structure given a set of nodal data.
  void mComputeMatrixStructure(const DataBase<Dimension>& dataBase,
                               const State<Dimension>& state) const;
  
  //! Update the Laplacian matrix given a set of nodal data (but don't alter its
  //! non-zero structure).
  void mUpdateLaplacianMatrix(const DataBase<Dimension>& dataBase,
                              const State<Dimension>& state) const;
  
  //! Create a right-hand-side vector for the linear system representing the 
  //! Poisson equation.
  PyObject* mCreateRHS(const State<Dimension>& state) const;
  
  //! Compute the gravitational potential given a set of nodal data.
  void mComputeGravitationalPotential(
          const DataBase<Dimension>& dataBase,
          const State<Dimension>& state) const;
  
  // Default constructor -- disabled.
  SPHGravity();

  // Copy constructor -- disabled.
  SPHGravity(const SPHGravity&);

  // Assignment operator -- disabled.
  SPHGravity& operator=(const SPHGravity&);

}; // end class SPHGravity

} // end namespace Spheral

#endif
