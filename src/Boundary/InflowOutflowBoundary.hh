//---------------------------------Spheral++----------------------------------//
// InflowOutflowBoundary -- creates inflow ghost images, which become internal nodes
// as they cross the specified boundary plane.
//
// Created by JMO, Tue Oct 15 11:23:09 PDT 2019
//
// Modified by:
//----------------------------------------------------------------------------//
#ifndef __Spheral_InflowOutflowBoundary__
#define __Spheral_InflowOutflowBoundary__

#include "Boundary.hh"
#include "Physics/Physics.hh"
#include "Geometry/GeomPlane.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/StateBase.hh" // For constructing Field keys.

namespace Spheral {

// Forward declarations.
template<typename Dimension> class NodeList;
template<typename Dimension> class FieldBase;
template<typename Dimension, typename DataType> class Field;
template<typename Dimension, typename DataType> class FieldList;
template<typename Dimension> class DataBase;

template<typename Dimension>
class InflowOutflowBoundary: public Boundary<Dimension>, public Physics<Dimension> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  typedef typename StateBase<Dimension>::KeyType KeyType;
  typedef typename Physics<Dimension>::TimeStepType TimeStepType;

  // Constructors and destructors.
  InflowOutflowBoundary(DataBase<Dimension>& dataBase,
                        const GeomPlane<Dimension>& plane,
                        const bool empty);
  virtual ~InflowOutflowBoundary() = default;

  //**********************************************************************
  // Boundary condition methods:
  // Use the given NodeList's neighbor object to select the ghost nodes.
  virtual void setGhostNodes(NodeList<Dimension>& nodeList) override;

  // For the computed set of ghost nodes, set the positions and H's.
  virtual void updateGhostNodes(NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the ghost node values in the given Field.
  virtual void applyGhostBoundary(FieldBase<Dimension>& field) const override;

  // Find any internal nodes that are in violation of this Boundary.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override;

  // For the computed set of nodes in violation of the boundary, bring them
  // back into compliance (for the positions and H's.)
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override;

  // This boundary does not cull ghosts, but others might have.
  virtual void cullGhostNodes(const FieldList<Dimension, size_t>& flagSet,
                              FieldList<Dimension, size_t>& old2newIndexMap,
                              std::vector<size_t>& numNodesRemoved) override;

  // After physics have been initialized we take a snapshot of the node state.
  virtual void initializeProblemStartup(const bool final) override;

  // We need to not cull ghost nodes, since they might need to cross the boundary
  // and become new inflow nodes.
  virtual bool allowGhostCulling() const override { return false; }
  //**********************************************************************

  //**********************************************************************
  // Physics methods:
  virtual void evaluateDerivatives(const Scalar time,
                                   const Scalar dt,
                                   const DataBase<Dimension>& dataBase,
                                   const State<Dimension>& state,
                                   StateDerivatives<Dimension>& derivatives) const override;

  // Vote on a time step.
  virtual TimeStepType dt(const DataBase<Dimension>& dataBase, 
                          const State<Dimension>& state,
                          const StateDerivatives<Dimension>& derivs,
                          const Scalar currentTime) const override;

  // Register the state you want carried around (and potentially evolved), as
  // well as the policies for such evolution.
  virtual void registerState(DataBase<Dimension>& dataBase,
                             State<Dimension>& state) override;

  // Register the derivatives/change fields for updating state.
  virtual void registerDerivatives(DataBase<Dimension>& dataBase,
                                   StateDerivatives<Dimension>& derivs) override;

  // Packages might want a hook to do some post-step finalizations.
  // Really we should rename this post-step finalize.
  virtual void finalize(const Scalar time, 
                        const Scalar dt,
                        DataBase<Dimension>& dataBase, 
                        State<Dimension>& state,
                        StateDerivatives<Dimension>& derivs) override;
  //**********************************************************************

  // Accessor methods.
  Scalar dtmin() const;
  const DataBase<Dimension>& dataBase() const;
  const GeomPlane<Dimension>& plane() const;
  size_t numInflowNodes(const NodeList<Dimension>& nodeList) const;
  
  void inflowRadius(const Scalar x);
  Scalar inflowRadius() const;

  // Get the stored data for generating ghost nodes.
  template<typename DataType> std::vector<DataType> storedValues(const KeyType key, const DataType& dummy);
  template<typename DataType> std::vector<DataType> storedValues(const Field<Dimension, DataType>& field);
  std::vector<std::string> storedKeys() const;

  // Set new values for the ghost nodes.
  template<typename DataType> void setStoredValues(const KeyType key, const std::vector<DataType>& values);
  template<typename DataType> void setStoredValues(const Field<Dimension, DataType>& field, const std::vector<DataType>& values);

  // Set new (constant) values for the ghost nodes.
  template<typename DataType> void setStoredValues(const KeyType key, const DataType& value);
  template<typename DataType> void setStoredValues(const Field<Dimension, DataType>& field, const DataType& value);

  // Clear out stored values.
  void clearStoredValues();

  //****************************************************************************
  // Methods required for restarting.
  virtual std::string label() const override;
  virtual void dumpState(FileIO& file, const std::string& pathName) const;
  virtual void restoreState(const FileIO& file, const std::string& pathName);
  //****************************************************************************

  // Prevent the Boundary virtual methods from being hidden
  using Boundary<Dimension>::applyGhostBoundary;
  using Boundary<Dimension>::enforceBoundary;
  using Physics<Dimension>::initializeProblemStartup;

private:
  //--------------------------- Private Interface ---------------------------//
  Scalar mInflowRadius;   // radius to clip inflow

  DataBase<Dimension>& mDataBase;
  GeomPlane<Dimension> mPlane;
  size_t mBoundaryCount;
  Scalar mDT;
  bool mActive, mEmpty;
  std::map<std::string, size_t> mNumInflowNodes;
  std::map<std::string, Scalar> mXmin;

  typedef std::map<KeyType, std::vector<char>> StorageType;
  StorageType mBufferedValues;

  // The restart registration.
  RestartRegistrationType mRestart;

  // No default or copy constructors.
  InflowOutflowBoundary();
  InflowOutflowBoundary(InflowOutflowBoundary&);
};

}

#include "InflowOutflowBoundaryInline.hh"

#endif
