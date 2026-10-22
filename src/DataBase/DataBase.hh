//---------------------------------Spheral++----------------------------------//
// DataBase -- The central point to store NodeLists and FieldLists, to serve the 
//             Spheral++ physics modules.
//
// Created by JMO, Sun Feb  6 13:44:49 PST 2000
//----------------------------------------------------------------------------//
#ifndef __Spheral_DataBase__
#define __Spheral_DataBase__

#include "DataBase/DataBaseView.hh"
#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "NodeList/SolidNodeList.hh"
#include "NodeList/DEMNodeList.hh"
#include "Neighbor/ConnectivityMap.hh"

#include <vector>
#include <span>

namespace Spheral {
  template<typename Dimension, typename DataType> class Field;
  template<typename Dimension, typename DataType> class FieldList;
  template<typename Dimension> class TableKernel;
}

namespace Spheral {

template<typename Dimension>
class DataBase: public DataBaseView<Dimension> {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;

  template<typename T> using ValueContainer = std::vector<T>;

#ifdef SPHERAL_UNIFIED_MEMORY
  template<typename T> using ViewSpan = std::span<T>;
#else
  template<typename T> using ViewSpan = chai::ManagedArray<T>;
#endif

  using NodeListIterator = typename ValueContainer<NodeList<Dimension>*>::iterator;
  using ConstNodeListIterator = typename ValueContainer<NodeList<Dimension>*>::const_iterator;

  using FluidNodeListIterator = typename ValueContainer<FluidNodeList<Dimension>*>::iterator;
  using ConstFluidNodeListIterator = typename ValueContainer<FluidNodeList<Dimension>*>::const_iterator;

  using SolidNodeListIterator = typename ValueContainer<SolidNodeList<Dimension>*>::iterator;
  using ConstSolidNodeListIterator = typename ValueContainer<SolidNodeList<Dimension>*>::const_iterator;

  using DEMNodeListIterator = typename ValueContainer<DEMNodeList<Dimension>*>::iterator;
  using ConstDEMNodeListIterator = typename ValueContainer<DEMNodeList<Dimension>*>::const_iterator;

  using ConnectivityMapType = ConnectivityMap<Dimension>;
  using ConnectivityMapPtr = std::shared_ptr<ConnectivityMapType>;
  
  using ViewType = DataBaseView<Dimension>;

  // Constructors, destructor
  DataBase();
  DataBase(const DataBase& rhs) = default;
  DataBase(DataBase&& rhs) = default;
  virtual ~DataBase() = default;

  // Assignment.
  DataBase& operator=(const DataBase& rhs) = default;
  DataBase& operator=(DataBase&& rhs) = default;

  // Global numbers of nodes
  size_t globalNumInternalNodes() const;
  size_t globalNumGhostNodes() const;
  size_t globalNumNodes() const;

  // Global numbers of fluid nodes
  size_t globalNumFluidInternalNodes() const;
  size_t globalNumFluidGhostNodes() const;
  size_t globalNumFluidNodes() const;
   
  // Global numbers of solid nodes
  size_t globalNumSolidInternalNodes() const;
  size_t globalNumSolidGhostNodes() const;
  size_t globalNumSolidNodes() const;
   
  // Global numbers of DEM nodes
  size_t globalNumDEMInternalNodes() const;
  size_t globalNumDEMGhostNodes() const;
  size_t globalNumDEMNodes() const;
   
  // Access the various NodeLists
  const ValueContainer<NodeList<Dimension>*>& nodeListPtrs()           const { return mNodeListPtrs; }
  const ValueContainer<FluidNodeList<Dimension>*>& fluidNodeListPtrs() const { return mFluidNodeListPtrs; }
  const ValueContainer<SolidNodeList<Dimension>*>& solidNodeListPtrs() const { return mSolidNodeListPtrs; }
  const ValueContainer<DEMNodeList<Dimension>*>& DEMNodeListPtrs()     const { return mDEMNodeListPtrs; }

  // Provide normal iterator methods over the DataBase NodeLists.
  NodeListIterator nodeListBegin()                                           { return mNodeListPtrs.begin(); }
  NodeListIterator nodeListEnd()                                             { return mNodeListPtrs.end(); }

  ConstNodeListIterator nodeListBegin()                                const { return mNodeListPtrs.begin(); }
  ConstNodeListIterator nodeListEnd()                                  const { return mNodeListPtrs.end(); }

  FluidNodeListIterator fluidNodeListBegin()                                 { return mFluidNodeListPtrs.begin(); }
  FluidNodeListIterator fluidNodeListEnd()                                   { return mFluidNodeListPtrs.end(); }

  ConstFluidNodeListIterator fluidNodeListBegin()                      const { return mFluidNodeListPtrs.begin(); }
  ConstFluidNodeListIterator fluidNodeListEnd()                        const { return mFluidNodeListPtrs.end(); }

  NodeListIterator fluidNodeListAsNodeListBegin()                            { return mFluidNodeListAsNodeListPtrs.begin(); }
  NodeListIterator fluidNodeListAsNodeListEnd()                              { return mFluidNodeListAsNodeListPtrs.end(); }

  ConstNodeListIterator fluidNodeListAsNodeListBegin()                 const { return mFluidNodeListAsNodeListPtrs.begin(); }
  ConstNodeListIterator fluidNodeListAsNodeListEnd()                   const { return mFluidNodeListAsNodeListPtrs.end(); }

  SolidNodeListIterator solidNodeListBegin()                                 { return mSolidNodeListPtrs.begin(); }
  SolidNodeListIterator solidNodeListEnd()                                   { return mSolidNodeListPtrs.end(); }

  ConstSolidNodeListIterator solidNodeListBegin()                      const { return mSolidNodeListPtrs.begin(); }
  ConstSolidNodeListIterator solidNodeListEnd()                        const { return mSolidNodeListPtrs.end(); }

  NodeListIterator solidNodeListAsNodeListBegin()                            { return mSolidNodeListAsNodeListPtrs.begin(); }
  NodeListIterator solidNodeListAsNodeListEnd()                              { return mSolidNodeListAsNodeListPtrs.end(); }

  ConstNodeListIterator solidNodeListAsNodeListBegin()                 const { return mSolidNodeListAsNodeListPtrs.begin(); }
  ConstNodeListIterator solidNodeListAsNodeListEnd()                   const { return mSolidNodeListAsNodeListPtrs.end(); }
  
  DEMNodeListIterator DEMNodeListBegin()                                     { return mDEMNodeListPtrs.begin(); }
  DEMNodeListIterator DEMNodeListEnd()                                       { return mDEMNodeListPtrs.end(); }

  ConstDEMNodeListIterator DEMNodeListBegin()                          const { return mDEMNodeListPtrs.begin(); }
  ConstDEMNodeListIterator DEMNodeListEnd()                            const { return mDEMNodeListPtrs.end(); }

  NodeListIterator DEMNodeListAsNodeListBegin()                              { return mDEMNodeListAsNodeListPtrs.begin(); }
  NodeListIterator DEMNodeListAsNodeListEnd()                                { return mDEMNodeListAsNodeListPtrs.end(); }

  ConstNodeListIterator DEMNodeListAsNodeListBegin()                   const { return mDEMNodeListAsNodeListPtrs.begin(); }
  ConstNodeListIterator DEMNodeListAsNodeListEnd()                     const { return mDEMNodeListAsNodeListPtrs.end(); }

  // Optimize all Neighbor objects for the current state.
  void reinitializeNeighbors() const;

  // Update the internal connectivity map.
  void updateConnectivityMap(const bool computeGhostConnectivity,
                             const bool computeOverlapConnectivity,
                             const bool computeIntersectionConnectivity) const;
  void patchConnectivityMap(const FieldList<Dimension, size_t>& flags,
                            const FieldList<Dimension, size_t>& old2new) const;

  // Get the connectivity map.
  const ConnectivityMapType& connectivityMap() const;
  ConnectivityMapPtr connectivityMapPtr() const;
  const ConnectivityMapType& connectivityMap(const bool computeGhostConnectivity,
                                             const bool computeOverlapConnectivity,
                                             const bool computeIntersectionConnectivity) const;
  ConnectivityMapPtr connectivityMapPtr(const bool computeGhostConnectivity,
                                        const bool computeOverlapConnectivity,
                                        const bool computeIntersectionConnectivity) const;

  // Methods to add, remove, and verify NodeLists.
  void appendNodeList(DEMNodeList<Dimension>& nodeList);
  void appendNodeList(SolidNodeList<Dimension>& nodeList);
  void appendNodeList(FluidNodeList<Dimension>& nodeList);
  void appendNodeList(NodeList<Dimension>& nodeList);

  void deleteNodeList(DEMNodeList<Dimension>& nodeList);
  void deleteNodeList(SolidNodeList<Dimension>& nodeList);
  void deleteNodeList(FluidNodeList<Dimension>& nodeList);
  void deleteNodeList(NodeList<Dimension>& nodeList);

  // Query methods which return "global" fields (FieldLists) for quantities
  // defined over all NodeLists.  These methods all build up FieldLists
  // referencing internally cached data in NodeLists, and so are fast.
  FieldList<Dimension, Scalar> globalMass() const;
  FieldList<Dimension, Scalar> globalMassRZ() const;
  FieldList<Dimension, Vector> globalPosition() const;
  FieldList<Dimension, Vector> globalVelocity() const;
  FieldList<Dimension, SymTensor> globalHfield() const;
  FieldList<Dimension, Scalar> globalWork() const;

  FieldList<Dimension, Scalar> fluidMass() const;
  FieldList<Dimension, Scalar> fluidMassRZ() const;
  FieldList<Dimension, Vector> fluidPosition() const;
  FieldList<Dimension, Vector> fluidVelocity() const;
  FieldList<Dimension, Scalar> fluidMassDensity() const;
  FieldList<Dimension, Scalar> fluidMassDensityRZ() const;
  FieldList<Dimension, Scalar> fluidSpecificThermalEnergy() const;
  FieldList<Dimension, SymTensor> fluidHfield() const;
  FieldList<Dimension, Scalar> fluidWork() const;

  FieldList<Dimension, Scalar> solidMass() const;
  FieldList<Dimension, Scalar> solidMassRZ() const;
  FieldList<Dimension, Vector> solidPosition() const;
  FieldList<Dimension, Vector> solidVelocity() const;
  FieldList<Dimension, Scalar> solidMassDensity() const;
  FieldList<Dimension, Scalar> solidMassDensityRZ() const;
  FieldList<Dimension, Scalar> solidSpecificThermalEnergy() const;
  FieldList<Dimension, SymTensor> solidHfield() const;
  FieldList<Dimension, Scalar> solidWork() const;
  FieldList<Dimension, SymTensor> solidDeviatoricStress() const;
  FieldList<Dimension, Scalar> solidPlasticStrain() const;
  FieldList<Dimension, Scalar> solidPlasticStrainRate() const;
  FieldList<Dimension, SymTensor> solidDamage() const;
  FieldList<Dimension, int> solidFragmentIDs() const;
  FieldList<Dimension, int> solidParticleTypes() const;

  FieldList<Dimension, Scalar> DEMMass() const;
  FieldList<Dimension, Vector> DEMPosition() const;
  FieldList<Dimension, Vector> DEMVelocity() const;
  FieldList<Dimension, SymTensor> DEMHfield() const;
  FieldList<Dimension, Scalar> DEMParticleRadius() const;
  FieldList<Dimension, int> DEMCompositeParticleIndex() const;
  FieldList<Dimension, size_t> DEMUniqueIndex() const;

  void setDEMUniqueIndices();

  // We can also return the node extent Fields stored in the Neighbor objects.
  FieldList<Dimension, Vector> globalNodeExtent() const;
  FieldList<Dimension, Vector> fluidNodeExtent() const;
  FieldList<Dimension, Vector> solidNodeExtent() const;
  FieldList<Dimension, Vector> DEMNodeExtent() const;

  // These functions return FieldLists with Fields that have to be calculated and
  // stored, so they are more expensive.
  // These FieldLists represent state that must be computed on the fly.
  void globalHinverse(FieldList<Dimension, SymTensor>& result) const;
  void fluidHinverse(FieldList<Dimension, SymTensor>& result) const;
  void fluidPressure(FieldList<Dimension, Scalar>& result) const;
  void fluidTemperature(FieldList<Dimension, Scalar>& result) const;
  void fluidSoundSpeed(FieldList<Dimension, Scalar>& result) const;
  void fluidVolume(FieldList<Dimension, Scalar>& result) const;
  void fluidGamma(FieldList<Dimension, Scalar>& result) const;
  void fluidEntropy(FieldList<Dimension, Scalar>& result) const;
  void fluidLinearMomentum(FieldList<Dimension, Vector>& result) const;
  void fluidTotalEnergy(FieldList<Dimension, Scalar>& result) const;
  void fluidSpecificHeat(const FieldList<Dimension, Scalar>& temperature,
                         FieldList<Dimension, Scalar>& result) const;

  // Collect the number of neighbors for each node from the ConnectivityMap.
  FieldList<Dimension, size_t> numNeighbors() const;

  // Return the number of internal nodes for the sets of NodeLists
  std::vector<size_t> numInternalNodesPerNodeList() const;
  std::vector<size_t> numInternalNodesPerFluidNodeList() const;
  std::vector<size_t> numInternalNodesPerSolidNodeList() const;
  std::vector<size_t> numInternalNodesPerDEMNodeList() const;

  //............................................................................
  // Create new FieldLists of size the number of NodeLists or FluidNodeLists.
  template<typename DataType>
  FieldList<Dimension, DataType> newGlobalFieldList(const DataType value,
                                                    const typename Field<Dimension, DataType>::FieldName name = "Unnamed Field") const;
  template<typename DataType>
  FieldList<Dimension, DataType> newFluidFieldList(const DataType value,
                                                   const typename Field<Dimension, DataType>::FieldName name = "Unnamed Field") const;
  template<typename DataType>
  FieldList<Dimension, DataType> newSolidFieldList(const DataType value,
                                                   const typename Field<Dimension, DataType>::FieldName name = "Unnamed Field") const;
  template<typename DataType>
  FieldList<Dimension, DataType> newDEMFieldList(const DataType value,
                                                   const typename Field<Dimension, DataType>::FieldName name = "Unnamed Field") const;

  // Resize a FieldList to the number of NodeLists or FluidNodeLists.
  // Optionally we can also set all elements in the FieldList to the specified value.
  // Note that if the FieldList is resized it is reconstructed from scratch, so all elements
  // will get the new value regardless of resetValues.
  template<typename DataType>
  void resizeGlobalFieldList(FieldList<Dimension, DataType>& fieldList,
                             const DataType value,
                             const typename Field<Dimension, DataType>::FieldName name = "Unnamed Field",
                             const bool resetValues = true) const;
  template<typename DataType>
  void resizeFluidFieldList(FieldList<Dimension, DataType>& fieldList,
                            const DataType value,
                            const typename Field<Dimension, DataType>::FieldName name = "Unnamed Field",
                            const bool resetValues = true) const;
  template<typename DataType>
  void resizeSolidFieldList(FieldList<Dimension, DataType>& fieldList,
                            const DataType value,
                            const typename Field<Dimension, DataType>::FieldName name = "Unnamed Field",
                            const bool resetValues = true) const;
  template<typename DataType>
  void resizeDEMFieldList(FieldList<Dimension, DataType>& fieldList,
                            const DataType value,
                            const typename Field<Dimension, DataType>::FieldName name = "Unnamed Field",
                            const bool resetValues = true) const;

  //............................................................................
  // Create vector<vector<>> versions of the FieldLists.
  template<typename DataType> std::vector<std::vector<DataType>> newGlobalArray(const DataType value) const;
  template<typename DataType> std::vector<std::vector<DataType>> newFluidArray(const DataType value) const;
  template<typename DataType> std::vector<std::vector<DataType>> newSolidArray(const DataType value) const;
  template<typename DataType> std::vector<std::vector<DataType>> newDEMArray(const DataType value) const;

  // Resize vector<vector<>> versions of the FieldLists.
  template<typename DataType> void resizeGlobalArray(std::vector<std::vector<DataType>>& array,
                                                     const DataType value,
                                                     const bool resetValues = true) const;
  template<typename DataType> void resizeFluidArray(std::vector<std::vector<DataType>>& array,
                                                    const DataType value,
                                                    const bool resetValues = true) const;
  template<typename DataType> void resizeSolidArray(std::vector<std::vector<DataType>>& array,
                                                    const DataType value,
                                                    const bool resetValues = true) const;
  template<typename DataType> void resizeDEMArray(std::vector<std::vector<DataType>>& array,
                                                    const DataType value,
                                                    const bool resetValues = true) const;
  //............................................................................

  // Compute coordinates bounding all nodes in the DataBase.
  void boundingBox(Vector& xmin, Vector& xmax,
                   const bool ghost = true) const;
  void boundingBox(Vector& xmin, Vector& xmax,
                   const FieldList<Dimension, int>& mask,
                   const bool ghost = true) const;

  // Return the local and global max sampling extents.
  void localSamplingBoundingVolume(Vector& centroid, double& radiusNodes, double& radiusSample,
                                   Vector& xminNodes, Vector& xmaxNodes,
                                   Vector& xminSample, Vector& xmaxSample) const;
  void globalSamplingBoundingVolume(Vector& centroid, double& radiusNodes, double& radiusSample,
                                    Vector& xminNodes, Vector& xmaxNodes,
                                    Vector& xminSample, Vector& xmaxSample) const;

  // Return the local min and max sampling extents for groupings of connected nodes.
  void localSamplingBoundingBoxes(std::vector<Vector>& xminima,
                                  std::vector<Vector>& xmaxima) const;
  void globalSamplingBoundingBoxes(std::vector<Vector>& xminima,
                                   std::vector<Vector>& xmaxima) const;

  // Provide a method to determine if the DataBase is in a minimally defined
  // valid state.
  bool valid() const;

  // Manage the View object
  ViewType view()  { return static_cast<DataBaseView<Dimension>&>(*this); }

private:
  //--------------------------- Private Interface ---------------------------//
  // NodeList pointer arrays
  ValueContainer<NodeList<Dimension>*>      mNodeListPtrs;
  ValueContainer<FluidNodeList<Dimension>*> mFluidNodeListPtrs;
  ValueContainer<NodeList<Dimension>*>      mFluidNodeListAsNodeListPtrs;
  ValueContainer<SolidNodeList<Dimension>*> mSolidNodeListPtrs;
  ValueContainer<NodeList<Dimension>*>      mSolidNodeListAsNodeListPtrs;
  ValueContainer<DEMNodeList<Dimension>*>   mDEMNodeListPtrs;
  ValueContainer<NodeList<Dimension>*>      mDEMNodeListAsNodeListPtrs;

  // Neighbor stuff
  mutable ConnectivityMapPtr mConnectivityMapPtr;

  // Storage for DataBaseView
#ifdef SPHERAL_UNIFIED_MEMORY  
  ValueContainer<NodeListView<Dimension>>      mNodeListViewStorage;
  ValueContainer<FluidNodeListView<Dimension>> mFluidNodeListViewStorage;
  ValueContainer<SolidNodeListView<Dimension>> mSolidNodeListViewStorage;
  ValueContainer<DEMNodeListView<Dimension>>   mDEMNodeListViewStorage;
#else
  ViewSpan<NodeListView<Dimension>>      mNodeListViewStorage;
  ViewSpan<FluidNodeListView<Dimension>> mFluidNodeListViewStorage;
  ViewSpan<SolidNodeListView<Dimension>> mSolidNodeListViewStorage;
  ViewSpan<DEMNodeListView<Dimension>>   mDEMNodeListViewStorage;
#endif

  // Keep the DataBaseView in sync
  void updateView();
};

}

#include "DataBaseInline.hh"

#endif
