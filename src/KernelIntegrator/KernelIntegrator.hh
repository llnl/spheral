//---------------------------------Spheral++----------------------------------//
// KernelIntegrator
//
// Performs integrals of kernels and coefficients
//----------------------------------------------------------------------------//
#ifndef __Spheral_KernelIntegrator_hh__
#define __Spheral_KernelIntegrator_hh__

#include <memory>
#include <vector>
#include <unordered_map>
#include "FlatConnectivity.hh"
#include "KernelIntegral.hh"
#include "IntegrationKernel.hh"
#include "DataBase/State.hh"
#include "Geometry/Dimension.hh"
#include "RK/RKCorrectionParams.hh"

namespace Spheral {

template<typename Dimension>
class KernelIntegrator {
public:
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using SymTensor = typename Dimension::SymTensor;
  using Tensor = typename Dimension::Tensor;
  using FacetedVolume = typename Dimension::FacetedVolume;
  using Facet = typename Dimension::Facet;
  using Subfacet = typename std::array<Vector, Dimension::nDim>;
  using ArrayDim = typename std::array<int, Dimension::nDim>;
  
  // Construtor
  KernelIntegrator(const int integrationOrder,
                   const std::shared_ptr<IntegrationKernel<Dimension>> kernel,
                   const DataBase<Dimension>& dataBase,
                   const FlatConnectivity<Dimension>& flatConnectivity);

  // Set the state, which should include Voronoi information
  virtual void setState(const double time,
                        const State<Dimension>& state);
  
  // Add an integral
  virtual void addIntegral(std::shared_ptr<KernelIntegralBase<Dimension>> integral);
  
  // Perform the integration
  virtual void performIntegration();

  // Get the flat connectivity
  virtual const FlatConnectivity<Dimension>& getFlatConnectivity() const { return mFlatConnectivity; }

  // Given coefficients, return values at centers
  template<typename DataType>
  void coefficientsToValue(const FieldList<Dimension, DataType>& coeffs,
                           FieldList<Dimension, DataType>& value) const;

  // Get the number of subcells integrated over
  int totalNumSubcells() const { return mTotalNumSubcells; }

  // Get the number of subfacets integrated over
  int totalNumSubfacets() const { return mTotalNumSubfacets; }

private:
  // Initialize the quadrature
  void initializeQuadrature();
  
  // // Get the subcells
  // // Top-level calls should have the subcells empty
  // void getSubcells(const FacetedVolume& cell,
  //                  std::vector<FacetedVolume>& subcells);
  
  // Get a quadrature
  void getQuadrature(const FacetedVolume& region,
                     std::vector<Scalar>& weights,
                     std::vector<Vector>& ordinates);

  // // Get the surface subcells
  // void getSubfacets(const Facet& facet,
  //                   std::vector<Subfacet>& subfacets);
  
  // Get surface indices
  void getSurfaceIndices(const ArrayDim& normal,
                         const std::vector<int>& indices,
                         std::vector<int>& surfaceIndex,
                         std::vector<int>& numSurfaces);
  
  // Get a surface quadrature
  void getSurfaceQuadrature(const Subfacet& subfacet,
                            std::vector<Scalar>& weights,
                            std::vector<Vector>& ordinates);

  // Input data
  const int mIntegrationOrder;
  const std::shared_ptr<IntegrationKernel<Dimension>> mKernel;
  const DataBase<Dimension>& mDataBase;
  const FlatConnectivity<Dimension>& mFlatConnectivity;

  // Constant for the kernel extent
  const double mHmult;
  
  // Base quadratures
  int mNumOrdinates;
  std::vector<Scalar> mBaseWeights;
  std::vector<Vector> mBaseOrdinates;

  // Surface quadratures
  int mNumSurfaceOrdinates;
  std::vector<Scalar> mBaseSurfaceWeights;
  std::vector<Vector> mBaseSurfaceOrdinates;
  
  // State
  double mTime;
  std::unique_ptr<State<Dimension>> mState;

  // Integrals that we want to add to
  std::vector<std::shared_ptr<KernelIntegralBase<Dimension>>> mIntegrals;

  // Scratch data, if we need it for evaluating coefficients
  mutable KernelIntegrationData<Dimension> mScratchData;

  // Keep track of the number of cells and faces we integrate over
  mutable int mTotalNumSubcells;
  mutable int mTotalNumSubfacets;
  
}; // end class KernelIntegrator

} // end namespace Spheral

#include "KernelIntegratorInline.hh"

#endif
