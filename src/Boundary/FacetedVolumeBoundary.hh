//---------------------------------Spheral++----------------------------------//
// FacetedVolumeBoundary -- Apply a FacetedVolume boundary condition to Spheral++
// Fields.
//
// Created by JMO, Sat Nov 30 21:12:34 PST 2019
//----------------------------------------------------------------------------//

#ifndef __Spheral_FacetedVolumeBoundary_hh__
#define __Spheral_FacetedVolumeBoundary_hh__

#include "Boundary.hh"

#include <vector>
#include <map>

namespace Spheral {

template<typename Dimension>
class FacetedVolumeBoundary: public Boundary<Dimension> {

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
  typedef typename FacetedVolume::Facet Facet;
  typedef GeomPlane<Dimension> Plane;

  // Constructors and destructors.
  FacetedVolumeBoundary(const FacetedVolume& poly,
                        const bool interiorBoundary,
                        const bool useGhosts);
  virtual ~FacetedVolumeBoundary();

  // Create any ghost nodes for the NodeList
  virtual void setGhostNodes(NodeList<Dimension>& nodeList) override;

  // For the computed set of ghost nodes, set the positions and H's.
  virtual void updateGhostNodes(NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the ghost values of given Field.
  virtual void applyGhostBoundary(Field<Dimension, int>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Scalar>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Vector>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, Tensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, SymTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FourthRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FifthRankTensor>& field) const override;
  virtual void applyGhostBoundary(Field<Dimension, FacetedVolume>& field) const override;

  // Find any internal nodes that are in violation of this Boundary.
  virtual void setViolationNodes(NodeList<Dimension>& nodeList) override;

  // For the computed set of nodes in violation of the boundary, bring them
  // back into compliance (for the positions and H's.)
  virtual void updateViolationNodes(NodeList<Dimension>& nodeList) override;

  // Apply the boundary condition to the violation node values in the given Field.
  virtual void enforceBoundary(Field<Dimension, int>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Scalar>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Vector>& field) const override;
  virtual void enforceBoundary(Field<Dimension, Tensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, SymTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, ThirdRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FourthRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FifthRankTensor>& field) const override;
  virtual void enforceBoundary(Field<Dimension, FacetedVolume>& field) const override;

  // Overridable hook for clearing out the boundary condition.
  virtual void reset(const DataBase<Dimension>& dataBase) override;

  // Use a set of flags to cull out inactive ghost nodes.
  virtual void cullGhostNodes(const FieldList<Dimension, size_t>& flagSet,
                              FieldList<Dimension, size_t>& old2newIndexMap,
                              std::vector<size_t>& numNodesRemoved) override;

  // Read access to internal data
  const FacetedVolume& polyVolume() const;
  bool interiorBoundary() const;
  bool useGhosts() const;
  const Tensor& reflectOperator(unsigned facetID) const;

  // Prevent the Boundary virtual methods from being hidden
  using Boundary<Dimension>::applyGhostBoundary;
  using Boundary<Dimension>::enforceBoundary;

private:
  //--------------------------- Private Interface ---------------------------//
  const FacetedVolume& mPoly;
  bool mInteriorBoundary, mUseGhosts;
  std::vector<Tensor> mReflectOperators;
  std::map<std::string, std::vector<std::vector<int>>> mFacetControlNodes;
  std::map<std::string, std::vector<std::pair<int,int>>> mFacetGhostNodes;
  std::map<std::string, std::vector<Tensor>> mViolationOperators;
};

}

#include "FacetedVolumeBoundaryInline.hh"

#endif
