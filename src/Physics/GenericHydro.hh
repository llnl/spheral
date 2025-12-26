//---------------------------------Spheral++----------------------------------//
// GenericHydro -- The base class for all Spheral++ hydro implementations.
//
// Created by JMO, Sat May 20 22:50:20 PDT 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_GenericHydro__
#define __Spheral_GenericHydro__

#include "Physics.hh"
#include "Geometry/Dimension.hh"
#include "SmoothingScale/SmoothingScaleBase.hh"

namespace Spheral {

template<typename Dimension> class ArtificialViscosity;
template<typename Dimension> class DataBase;

// Many hydro algorithms have these sorts of choices for the mass density and H.
enum class MassDensityType {
  SumDensity = 0,
  RigorousSumDensity = 1,
  HybridSumDensity = 2,
  IntegrateDensity = 3,
  VoronoiCellDensity = 4,
  SumVoronoiCellDensity = 5,
  CorrectedSumDensity = 6,
};

template<typename Dimension>
class GenericHydro: public Physics<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  using TimeStepType = typename Physics<Dimension>::TimeStepType;
  using ResidualType = typename Physics<Dimension>::ResidualType;

  // Constructors.
  GenericHydro(ArtificialViscosity<Dimension>& Q,
               const double cfl,
               const bool useVelocityMagnitudeForDt);

  // Destructor.
  virtual ~GenericHydro();

  // Forbidden methods.
  GenericHydro() = delete;
  GenericHydro(const GenericHydro&) = delete;
  GenericHydro& operator=(const GenericHydro&) = delete;

  // We require all Physics packages to provide a method returning their vote
  // for the next time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase,
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // Provide an appropriate vote for an implicit time step limit.
  virtual TimeStepType dtImplicit(const DataBase<Dimension>& dataBase,
                                  const State<Dimension>& state,
                                  const StateDerivatives<Dimension>& derivs,
                                  const Scalar currentTime) const override;

  // Return the maximum state change we care about for checking for convergence in the implicit integration methods.
  virtual ResidualType maxResidual(const DataBase<Dimension>& dataBase, 
                                   const State<Dimension>& state1,
                                   const State<Dimension>& state0,
                                   const Scalar tol) const override;

  // Allow access to the artificial viscosity.
  ArtificialViscosity<Dimension>& artificialViscosity() const { return mArtificialViscosity; }

  // Also allow access to the CFL timestep safety criteria.
  Scalar cfl()                              const { return mCFL; }
  void cfl(Scalar x)                              { mCFL = x; }

  // Attribute to determine if the absolute magnitude of the velocity should
  // be used in determining the timestep.
  bool useVelocityMagnitudeForDt()          const { return mUseVelocityMagnitudeForDt; }
  void useVelocityMagnitudeForDt(bool x)          { mUseVelocityMagnitudeForDt = x; }

  // Return the cumulative neighboring statistics.
  int minMasterNeighbor()                   const { return mMinMasterNeighbor; }
  int maxMasterNeighbor()                   const { return mMaxMasterNeighbor; }
  double averageMasterNeighbor()            const { return double(mSumMasterNeighbor)/(mNormMasterNeighbor + 1.0e-20); }

  int minCoarseNeighbor()                   const { return mMinCoarseNeighbor; }
  int maxCoarseNeighbor()                   const { return mMaxCoarseNeighbor; }
  double averageCoarseNeighbor()            const { return double(mSumCoarseNeighbor)/(mNormCoarseNeighbor + 1.0e-20); }

  int minRefineNeighbor()                   const { return mMinRefineNeighbor; }
  int maxRefineNeighbor()                   const { return mMaxRefineNeighbor; }
  double averageRefineNeighbor()            const { return double(mSumRefineNeighbor)/(mNormRefineNeighbor + 1.0e-20); }

  int minActualNeighbor()                   const { return mMinActualNeighbor; }
  int maxActualNeighbor()                   const { return mMaxActualNeighbor; }
  double averageActualNeighbor()            const { return double(mSumActualNeighbor)/(mNormActualNeighbor + 1.0e-20); }

  // Stored attributes about the last timestep chosen
  size_t DTrank()                           const { return mDTrank; }
  size_t DTNodeList()                       const { return mDTNodeList; }
  size_t DTnode()                           const { return mDTnode; }
  std::string DTreason()                    const { return mDTreason; }

  // Same things for residuals used in implicit time advancement
  size_t maxResidualRank()                  const { return mMaxResidualRank; }
  size_t maxResidualNodeList()              const { return mMaxResidualNodeList; }
  size_t maxResidualNode()                  const { return mMaxResidualNode; }
  std::string maxResidualReason()           const { return mMaxResidualReason; }

protected:
  //-------------------------- Protected Interface --------------------------//
  void updateMasterNeighborStats(int numMaster) const;
  void updateCoarseNeighborStats(int numNeighbor) const;
  void updateRefineNeighborStats(int numNeighbor) const;
  void updateActualNeighborStats(int numNeighbor) const;

private:
  //--------------------------- Private Interface ---------------------------//
  ArtificialViscosity<Dimension>& mArtificialViscosity;
  Scalar mCFL;
  bool mUseVelocityMagnitudeForDt;

  mutable int mMinMasterNeighbor, mMaxMasterNeighbor, mSumMasterNeighbor;
  mutable int mMinCoarseNeighbor, mMaxCoarseNeighbor, mSumCoarseNeighbor;
  mutable int mMinRefineNeighbor, mMaxRefineNeighbor, mSumRefineNeighbor;
  mutable int mMinActualNeighbor, mMaxActualNeighbor, mSumActualNeighbor;
  mutable int mNormMasterNeighbor;
  mutable int mNormCoarseNeighbor;
  mutable int mNormRefineNeighbor;
  mutable int mNormActualNeighbor;
  mutable size_t mDTrank, mDTNodeList, mDTnode;
  mutable std::string mDTreason;
  mutable size_t mMaxResidualRank, mMaxResidualNodeList, mMaxResidualNode;
  mutable std::string mMaxResidualReason;
};

}

#endif
