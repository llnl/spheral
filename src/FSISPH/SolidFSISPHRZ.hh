//---------------------------------Spheral++----------------------------------//
// SolidFSISPHRZ
//
// RZ version (2D cylindrical geometry) of SolidFSISPH
//
// This RZ version is a naive area-weighting implementation,
// similar to RZ versions of SPH and CRKSPH.
//----------------------------------------------------------------------------//
#ifndef __Spheral_SolidFSISPHRZ_hh__
#define __Spheral_SolidFSISPHRZ_hh__

#include "FSISPH/SolidFSISPH.hh"
#include "Geometry/Dimension.hh"

#include <string>
#include <vector>
#include <memory>

namespace Spheral {

class SolidFSISPHRZ: public SolidFSISPH<Dim<2>> {

public:

  //--------------------------- Public Interface ---------------------------//
  using Dimension = Dim<2>;
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using ConstBoundaryIterator = typename Physics<Dimension>::ConstBoundaryIterator;
  using PairAccelerationsType = PairwiseField<Dimension, Vector, 1u>;
  using PairWorkType = PairwiseField<Dimension, Scalar, 2u>;

  // Constructors.
  SolidFSISPHRZ(DataBase<Dimension>& dataBase,
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
  SolidFSISPHRZ() = delete;
  SolidFSISPHRZ(const SolidFSISPHRZ&) = delete;
  SolidFSISPHRZ& operator=(const SolidFSISPHRZ&) = delete;

  virtual ~SolidFSISPHRZ() = default;

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
  
  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override { return "SolidFSISPHRZ"; }
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
 //****************************************************************************
};

}

#endif
