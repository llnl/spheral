//---------------------------------Spheral++----------------------------------//
// VolumeUpdate
//
// Base class for volume computation packages. Owns the volume FieldList and
// dispatches to different volume computation methods. VoronoiCells inherits
// from this for full Voronoi geometry computation.
//----------------------------------------------------------------------------//
#ifndef __Spheral_VolumeUpdate__
#define __Spheral_VolumeUpdate__

#include "VoronoiCells/VolumeType.hh"
#include "DataOutput/registerWithRestart.hh"
#include "Field/FieldList.hh"
#include "Physics/Physics.hh"

namespace Spheral {
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class DataBase;
template<typename Dimension> class TableKernel;

template<typename Dimension>
class VolumeUpdate : public Physics<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;
  using TimeStepType = typename std::pair<double, std::string>;

  // Constructor
  VolumeUpdate(const VolumeType volumeType,
               const TableKernel<Dimension>& W,
               const bool updateInStep,
               const bool updateInFinalize);

  // Destructor
  virtual ~VolumeUpdate();

  // Volume computation — dispatches on mVolumeType.
  virtual void computeVolume(const DataBase<Dimension>& dataBase,
                             State<Dimension>& state);

  // Physics interface
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override;

  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  virtual std::string label() const override { return "VolumeUpdate"; }

  virtual void applyGhostBoundaries(State<Dimension>& state,
                                    StateDerivatives<Dimension>& derivs) override;

  virtual void enforceBoundaries(State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  virtual void initializeProblemStartup(DataBase<Dimension>& dataBase) override;

  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;

  virtual void preStepInitialize(const DataBase<Dimension>& dataBase,
                                 State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) override;

  virtual bool postStateUpdate(const Scalar time,
                               const Scalar dt,
                               const DataBase<Dimension>& dataBase,
                               State<Dimension>& state,
                               StateDerivatives<Dimension>& derivs) override;

  virtual void finalize(const Scalar time,
                                const Scalar dt,
                                DataBase<Dimension>& dataBase,
                                State<Dimension>& state,
                                StateDerivatives<Dimension>& derivs) override;

  // Accessors
  VolumeType volumeType()     const { return mVolumeType; }
  bool       updateInStep()   const { return mUpdateInStep; }
  bool       updateInFinalize() const { return mUpdateInFinalize; }
  const FieldList<Dimension, Scalar>& volume() const { return mVolume; }
  const FieldList<Dimension, Scalar>& volume3d() const { return mVolume3d; }

  // Restart
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);

protected:
  //--------------------------- Protected Interface ---------------------------//
  FieldList<Dimension, Scalar> mVolume;
  FieldList<Dimension, Scalar> mVolume3d;
  VolumeType mVolumeType;
  bool mUpdateInStep;
  bool mUpdateInFinalize;
  const TableKernel<Dimension>& mKernel;

private:
  //--------------------------- Private Interface ---------------------------//
  RestartRegistrationType mRestart;

  // No default constructor, copying, or assignment.
  VolumeUpdate();
  VolumeUpdate(const VolumeUpdate&);
  VolumeUpdate& operator=(const VolumeUpdate&);
}; // end VolumeUpdate

} // end namespace Spheral



#endif
