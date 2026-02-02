//---------------------------------Spheral++----------------------------------//
// SolidFSISPH -- SolidSPHHydro modified to better handle 
//                multimaterial material problems with interfaces
//                and large density ratios. 
//
// J.M. Pearl 2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidFSISPH_hh__
#define __Spheral_SolidFSISPH_hh__

#include "Physics/GenericHydro.hh"
#include "Utilities/SpheralMessage.hh"

#include <string>
#include <vector>
#include <memory>

namespace Spheral {

enum class InterfaceMethod {
  HLLCInterface = 0,
  ModulusInterface = 1,
  NoInterface = 2,
};

enum class KernelAveragingMethod {
  NeverAverageKernels = 0,
  AlwaysAverageKernels = 1,
  AverageInterfaceKernels = 2,
};

enum class FSIMassDensityMethod {
  FSISumMassDensity = 0,
  PressureCorrectSumMassDensity = 1,
  HWeightedSumMassDensity = 2,
};

template<typename Dimension> class State;
template<typename Dimension> class StateDerivatives;
template<typename Dimension> class ArtificialViscosityHandle;
template<typename Dimension> class SlideSurface;
template<typename Dimension> class TableKernel;
template<typename Dimension> class DataBase;
template<typename Dimension, typename Value> class Field;
template<typename Dimension, typename Value> class FieldList;
template<typename Dimension, typename Value, size_t numElements> class PairwiseField;
class FileIO;

template<typename Dimension>
class SolidFSISPH: public GenericHydro<Dimension> {

public:

  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;

  using PairAccelerationsType = PairwiseField<Dimension, Vector, 1u>;
  using PairWorkType = PairwiseField<Dimension, Scalar, 2u>;

  // Constructors.
  SolidFSISPH(DataBase<Dimension>& dataBase,
              ArtificialViscosityHandle<Dimension>& Q,
              SlideSurface<Dimension>& slide,
              const TableKernel<Dimension>& W,
              const double cfl,
              const double surfaceForceCoefficient,
              const double densityStabilizationCoefficient,
              const double specificThermalEnergyDiffusionCoefficient,
              const double xsphCoefficient,
              const InterfaceMethod interfaceMethod,
              const KernelAveragingMethod kernelAveragingMethod,
              const std::vector<int> sumDensityNodeLists,
              const bool useVelocityMagnitudeForDt,
              const bool compatibleEnergyEvolution,
              const bool evolveTotalEnergy,
              const bool linearCorrectGradients,
              const bool decoupleDamagedMaterial,
              const double interfacePmin,
              const double interfaceNeighborAngleThreshold,
              const FSIMassDensityMethod densityUpdate,
              const double epsTensile,
              const double nTensile,
              const Vector& xmin,
              const Vector& xmax);

  // No default constructor, copying, or assignment.
  SolidFSISPH() = delete;
  SolidFSISPH(const SolidFSISPH&) = delete;
  SolidFSISPH& operator=(const SolidFSISPH&) = delete;

  virtual ~SolidFSISPH() = default;

  // A second optional method to be called on startup, after Physics::initializeProblemStartup has
  // been called.
  // One use for this hook is to fill in dependendent state using the State object, such as
  // temperature or pressure.
  virtual void initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                                    State<Dimension>& state,
                                                    StateDerivatives<Dimension>& derivs) override;

  virtual
  void registerState(DataBase<Dimension>& dataBase,
                     State<Dimension>& state) override;

  virtual
  void registerDerivatives(DataBase<Dimension>& dataBase,
                           StateDerivatives<Dimension>& derivs) override;

  virtual 
  void preStepInitialize(const DataBase<Dimension>& dataBase, 
                               State<Dimension>& state,
                               StateDerivatives<Dimension>& derivs) override;
  virtual
  void evaluateDerivatives(const Scalar time,
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase,
                           const State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivatives) const override;
  void firstDerivativesLoop(const Scalar time,
                            const Scalar dt,
                            const DataBase<Dimension>& dataBase,
                            const State<Dimension>& state,
                                  StateDerivatives<Dimension>& derivatives) const;
  template<typename QType>
  void secondDerivativesLoop(const Scalar time,
                             const Scalar dt,
                             const DataBase<Dimension>& dataBase,
                             const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives,
                             const QType& Q) const;

  virtual 
  void finalizeDerivatives(const Scalar time, 
                           const Scalar dt,
                           const DataBase<Dimension>& dataBase, 
                           const State<Dimension>& state,
                                 StateDerivatives<Dimension>& derivs) const override;

  virtual
  void applyGhostBoundaries(State<Dimension>& state,
                            StateDerivatives<Dimension>& derivs) override;

  virtual
  void enforceBoundaries(State<Dimension>& state,
                         StateDerivatives<Dimension>& derivs) override;

  void linearReconstruction(const typename Dimension::Vector& ri,
                            const typename Dimension::Vector& rj,
                            const typename Dimension::Scalar& yi,
                            const typename Dimension::Scalar& yj,
                            const typename Dimension::Vector& DyDxi,
                            const typename Dimension::Vector& DyDxj,
                                  typename Dimension::Scalar& ytildei,
                                  typename Dimension::Scalar& ytildej) const;


  const TableKernel<Dimension>& kernel() const;
  SlideSurface<Dimension>& slideSurface() const;

  FSIMassDensityMethod densityUpdate() const;
  void densityUpdate(FSIMassDensityMethod type);

  InterfaceMethod interfaceMethod() const;
  void interfaceMethod(InterfaceMethod method);

  KernelAveragingMethod kernelAveragingMethod() const;
  void kernelAveragingMethod(KernelAveragingMethod method);

  bool compatibleEnergyEvolution() const;
  void compatibleEnergyEvolution(bool val);

  bool evolveTotalEnergy() const;
  void evolveTotalEnergy(bool val);

  bool linearCorrectGradients() const;
  void linearCorrectGradients(bool val);

  bool applySelectSumDensity() const;
  void applySelectSumDensity(bool x);

  bool planeStrain()                                const { DeprecationWarning("FSISPH WARNING: planeStrain is deprecated"); return false; }
  void planeStrain(bool val)                              { DeprecationWarning("FSISPH WARNING: planeStrain is deprecated"); }
  
  bool decoupleDamagedMaterial() const;
  void decoupleDamagedMaterial(bool val);

  std::vector<int> sumDensityNodeLists() const;
  void sumDensityNodeLists(std::vector<int> x);

  double surfaceForceCoefficient() const;
  void surfaceForceCoefficient(double x);

  double densityStabilizationCoefficient() const;
  void densityStabilizationCoefficient(double x);

  double specificThermalEnergyDiffusionCoefficient() const;
  void specificThermalEnergyDiffusionCoefficient(double x);

  double interfacePmin() const;
  void interfacePmin(double val);

  double interfaceNeighborAngleThreshold() const;
  void interfaceNeighborAngleThreshold(double val);

  double xsphCoefficient() const;
  void xsphCoefficient(double x);

  Scalar epsilonTensile() const;
  void epsilonTensile(Scalar val);

  Scalar nTensile() const;
  void nTensile(Scalar val);

  const Vector& xmin() const;
  const Vector& xmax() const;
  void xmin(const Vector& x);
  void xmax(const Vector& x);

  const PairAccelerationsType&           pairAccelerations() const;
  const PairWorkType&                    pairDepsDt() const;

  const FieldList<Dimension, int>&       timeStepMask() const;
  const FieldList<Dimension, Scalar>&    pressure() const;
  const FieldList<Dimension, Scalar>&    damagedPressure() const;
  const FieldList<Dimension, Scalar>&    soundSpeed() const;
  const FieldList<Dimension, Scalar>&    bulkModulus() const;
  const FieldList<Dimension, Scalar>&    shearModulus() const;
  const FieldList<Dimension, Scalar>&    yieldStrength() const;
  const FieldList<Dimension, Scalar>&    plasticStrain0() const;
  //const FieldList<Dimension, Scalar>&    inverseEquivalentDeviatoricStress() const;
  const FieldList<Dimension, Scalar>&    volume() const;
  const FieldList<Dimension, Vector>&    DxDt() const;
  const FieldList<Dimension, Vector>&    XSPHDeltaV() const;
  const FieldList<Dimension, Scalar>&    XSPHWeightSum() const;
  const FieldList<Dimension, Vector>&    DvDt() const;
  const FieldList<Dimension, Scalar>&    DmassDensityDt() const;
  const FieldList<Dimension, Scalar>&    DspecificThermalEnergyDt() const;
  const FieldList<Dimension, SymTensor>& DdeviatoricStressDt() const;
  const FieldList<Dimension, Vector>&    DPDx() const;
  const FieldList<Dimension, Vector>&    DepsDx() const;
  const FieldList<Dimension, Tensor>&    DvDx() const;
  const FieldList<Dimension, Tensor>&    internalDvDx() const;
  const FieldList<Dimension, Tensor>&    M() const;
  const FieldList<Dimension, Tensor>&    localM() const;
  const FieldList<Dimension, Scalar>&    normalization() const;

  const FieldList<Dimension, int>& interfaceFlags() const;
  const FieldList<Dimension, Vector>& interfaceAreaVectors() const;
  const FieldList<Dimension, Vector>& interfaceNormals() const;
  const FieldList<Dimension, Scalar>& interfaceSmoothness() const;

  const FieldList<Dimension, int>& newInterfaceFlags() const;
  const FieldList<Dimension, Vector>& newInterfaceAreaVectors() const;
  const FieldList<Dimension, Vector>& newInterfaceNormals() const;
  const FieldList<Dimension, Scalar>& interfaceSmoothnessNormalization() const;
  const FieldList<Dimension, Scalar>& interfaceFraction() const;
  const FieldList<Dimension, Scalar>& newInterfaceSmoothness() const;
  const FieldList<Dimension, Scalar>& interfaceAngles() const;
  
  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SolidFSISPH"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
 //****************************************************************************

private:
  const TableKernel<Dimension>& mKernel;
  SlideSurface<Dimension>& mSlideSurface;

  FSIMassDensityMethod mDensityUpdate;
  InterfaceMethod mInterfaceMethod;                   // switch for material interface method
  KernelAveragingMethod mKernelAveragingMethod;       // how do we handle our kernels?

  bool mCompatibleEnergyEvolution;
  bool mEvolveTotalEnergy; 
  bool mLinearCorrectGradients;
  bool mDecoupleDamagedMaterial;
  bool mApplySelectDensitySum;                        // switch for density sum
  std::vector<int> mSumDensityNodeLists;              // turn on density sum subset of nodeLists

  double mSurfaceForceCoefficient;                    // Monaghan 2013 force increase @ interface
  double mDensityStabilizationCoefficient;            // adjusts DvDx to stabilize rho
  double mSpecificThermalEnergyDiffusionCoefficient;  // controls diffusion of eps
  double mXSPHCoefficient;                            // controls amount of xsph-ing
  double mInterfacePmin;                              // min pressure across interfaces (similar to eos)
  double mInterfaceNeighborAngleThreshold;            // opening angle used to id surf/interface nodes

  Scalar mEpsTensile;
  Scalar mnTensile;

  Vector mxmin;
  Vector mxmax;

  std::unique_ptr<PairAccelerationsType> mPairAccelerationsPtr;    // store pairwise contribution to DvDt for compatible
  std::unique_ptr<PairWorkType>          mPairDepsDtPtr;           // store pairwise contribution to DepsDt for compatible

  FieldList<Dimension, int>       mTimeStepMask;
  FieldList<Dimension, Scalar>    mPressure;
  FieldList<Dimension, Scalar>    mDamagedPressure;
  FieldList<Dimension, Scalar>    mSoundSpeed;
  FieldList<Dimension, Scalar>    mBulkModulus;
  FieldList<Dimension, Scalar>    mShearModulus;
  FieldList<Dimension, Scalar>    mYieldStrength;
  FieldList<Dimension, Scalar>    mPlasticStrain0;
  //FieldList<Dimension, Scalar>    mInverseEquivalentDeviatoricStress;
  FieldList<Dimension, Scalar>    mVolume;
  FieldList<Dimension, Vector>    mDxDt;
  FieldList<Dimension, Vector>    mXSPHDeltaV;
  FieldList<Dimension, Scalar>    mXSPHWeightSum;
  FieldList<Dimension, Vector>    mDvDt;
  FieldList<Dimension, Scalar>    mDmassDensityDt;
  FieldList<Dimension, Scalar>    mDspecificThermalEnergyDt;
  FieldList<Dimension, SymTensor> mDdeviatoricStressDt;
  FieldList<Dimension, Vector>    mDPDx;
  FieldList<Dimension, Vector>    mDepsDx;
  FieldList<Dimension, Tensor>    mDvDx;
  FieldList<Dimension, Tensor>    mInternalDvDx;
  FieldList<Dimension, Tensor>    mM;
  FieldList<Dimension, Tensor>    mLocalM;
  FieldList<Dimension, Scalar>    mNormalization;

  FieldList<Dimension, int> mInterfaceFlags;                  // flags indicating interface type
  FieldList<Dimension, Vector> mInterfaceAreaVectors;         // interface area vectors that can be used for BCs
  FieldList<Dimension, Vector> mInterfaceNormals;             // surface normals between nodelists     
  FieldList<Dimension, Scalar> mInterfaceSmoothness;          // smoothness metric (0-1) 
  
  FieldList<Dimension, int> mNewInterfaceFlags;                   // flags indicating interface type
  FieldList<Dimension, Vector> mNewInterfaceAreaVectors;          // interface area vectors that can be used for BCs next time step
  FieldList<Dimension, Vector> mNewInterfaceNormals;              // surface normals between nodelists next time step    
  FieldList<Dimension, Scalar> mInterfaceSmoothnessNormalization; // normalization for a our smoothness metric
  FieldList<Dimension, Scalar> mInterfaceFraction;                // normalization for same material nodes
  FieldList<Dimension, Scalar> mNewInterfaceSmoothness;           // smoothness metric (0-1) next time step 
  FieldList<Dimension, Scalar> mInterfaceAngles;                  // check the angle for free-surface master nodes (type 2 -> type 3)

protected:
  //--------------------------- Protected Interface ---------------------------//
  // The restart registration.
  RestartRegistrationType mRestart;
};

}


#include "SolidFSISPHInline.hh"

#endif
