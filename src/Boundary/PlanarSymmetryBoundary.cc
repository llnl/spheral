#include "PlanarSymmetryBoundary.hh"
#include "NodeList/FluidNodeList.hh"

#include "Utilities/DBC.hh"

namespace Spheral {

using std::vector;

//------------------------------------------------------------------------------
PlanarSymmetryBoundary::
PlanarSymmetryBoundary(TableKernel<Dim<3> >* kernel):
   Boundary<Dim<3> >(),
   mKernel(kernel),
   mFirstGhostIndex(0),
   mL(0.0),
   mdz(0.0) {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
PlanarSymmetryBoundary::
~PlanarSymmetryBoundary() {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
setGhostNodes(NodeList<Dim<3> >& nodeList) {
   // Add this NodeList, ecreating space for control & ghost nodes.
   addNodeList(nodeList);

   // Define control nodes (and here, every internal node is a control node).
   Boundary<Dim<3> >::BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
   vector<int>& controlNodes = boundaryNodes.controlNodes;

   // We want to make walls of ghost nodes around the existing planar
   // internal nodes.  We will duplicate the XY plane along the Z axis over 
   // a distance L, which we define as the z-extent of the kernel for this 
   // boundary condition.

   // Now set up the ghost nodes for this boundary condition.
   // We create a bevy of nodes around each internal node.
   int numberOfGhostNodesToAdd = 0;
   const auto numInternalNodes = nodeList.numInternalNodes();
   for (auto nodei = 0u; nodei < numInternalNodes; ++nodei)
   {
      // Add this internal node as a control node.
      controlNodes.push_back(static_cast<int>(nodei));

      // Find the position of the internal node.
      Dim<3>::Vector xInternal = nodeList.positions()(nodei);

      // Now L is the total z-extent of the kernel, hz.  
      // Let M = H*H.  Then hz = sqrt(cofactor(M, 3, 3))/ det H.
      Dim<3>::SymTensor H = nodeList.Hfield()(nodei);
      double detH = H.Determinant();
      double M11 = H(0,0)*H(0,0) + H(0,1)*H(0,1) + H(0,2)*H(0,2);
      double M21 = H(0,0)*H(0,1) + H(0,1)*H(1,1) + H(0,2)*H(1,2);
//      double M31 = H(0,0)*H(0,2) + H(0,1)*H(1,2) + H(0,2)*H(2,2);
      double M22 = H(0,1)*H(0,1) + H(1,1)*H(1,1) + H(1,2)*H(1,2);
//      double M32 = H(0,1)*H(0,2) + H(1,1)*H(1,2) + H(1,2)*H(2,2);
//      double M33 = H(0,2)*H(0,2) + H(1,2)*H(1,2) + H(2,2)*H(2,2);
      double Mcof33 = M11*M22 - M21*M21;
      mL = Mcof33 / detH;

      // Now L = Nz*dz, where Nz is the number of neighbors in the z-extent of 
      // the kernel.  Thus, Nz = (# of neighbors per H) * (# of H per kernel z-extent).
      // Below, we add 1 to Nz to make sure that we're creating enough neighbors.
      // Note that dz = L / Nz.
      int Nz = static_cast<int>(nodeList.nodesPerSmoothingScale() * mKernel->kernelExtent()) + 1;
      mdz = mL / Nz;

      // We'll progress along the z direction, duplicating the XY plane.
      for (int i = 0; i < Nz; ++i) 
      {
         // Compute our z position.
         double z = i*mdz - (mL / 2.0);

         // Add the number of nodes in our computational XY plane.
         numberOfGhostNodesToAdd += numInternalNodes;
      } // end for
   } // end for

   // Allocate space for the new ghost nodes we're going to create.
   int currentNumGhostNodes = nodeList.numGhostNodes();
   nodeList.numGhostNodes(currentNumGhostNodes + numberOfGhostNodesToAdd);

   // Take note of the index of the first ghost node (and the last).
   mFirstGhostIndex = static_cast<size_t>(currentNumGhostNodes);

   // Let the neighbor object know that we've added ghost/control nodes.
   // FIXME: ????
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
updateGhostNodes(NodeList<Dim<3> >& nodeList) {
   // We are going to update all the positions of our ghost nodes.  We 
   // do this by iterating over the internal nodes and updating the 
   // ghost nodes that correspond to these internal nodes.

   // We must do this update only on the ghost nodes that pertain to this
   // particular boundary condition.
   const auto numInternalNodes = nodeList.numInternalNodes();
   auto ghostNodei = numInternalNodes + mFirstGhostIndex;

   // Now iterate over the ghost nodes of this boundary condition (there's one 
   // ghost node for every internal node) and move them along with the 
   // internal nodes.
   for (auto nodei = 0u; nodei < numInternalNodes; ++nodei)
   {
      // Find the position of the internal node.
      Dim<3>::Vector xInternal = nodeList.positions()(nodei);

      int Nz = static_cast<int>(mL / mdz) + 1;

      // We'll progress along the radial direction, sweeping out the azimuth.
      for (int i = 0; i < Nz; ++i) 
      {
         // Compute our z position.
         double z = i*mdz - (mL / 2.0);

         // Duplicate the positions of the XY plane.
         for (auto j = 0u; j < numInternalNodes; ++j)
         {
            // The x and y components of the ghost node are the same as that 
            // for the internal node.  The z component is given above.
            Dim<3>::Vector xGhost(xInternal(0), xInternal(1), z);
            nodeList.positions()(ghostNodei) = xGhost;
            ++ghostNodei;
         } // end for
      } // end for
   } // end for
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
setViolationNodes(NodeList<Dim<3> >& nodeList) {
   // Find the set of nodes in the given NodeList that violate this boundary
   // condition.  We flag all nodes that have moved out of the XY plane as 
   // violation nodes.
   Boundary<Dim<3> >::BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
   vector<int>& violationNodes = boundaryNodes.violationNodes;
   const auto numInternalNodes = nodeList.numInternalNodes();
   for (auto nodei = 0u; nodei < numInternalNodes; ++nodei)
   {
      violationNodes.push_back(static_cast<int>(nodei));
   } // end for
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
updateViolationNodes(NodeList<Dim<3> >& nodeList) {
   // Zero out the z components of every violation node.  We can do this
   // by looping through the violation nodes and just setting the positions 
   // of the nodes accordingly.
   Boundary<Dim<3> >::BoundaryNodes& boundaryNodes = accessBoundaryNodes(nodeList);
   vector<int>& violationNodes = boundaryNodes.violationNodes;
   Field<Dim<3>, Vector>& x = nodeList.positions();
   for (size_t i = 0; i < violationNodes.size(); ++i)
   {
      x[violationNodes[i]].z(0.0);
   } // end for
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Scalar>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Vector>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::Tensor>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
applyGhostBoundary(Field<Dim<3>, Dim<3>::SymTensor>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Scalar>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Vector>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
enforceBoundary(Field<Dim<3>, typename Dimension::Vector>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::Tensor>& field) const {
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
void
PlanarSymmetryBoundary::
enforceBoundary(Field<Dim<3>, Dim<3>::SymTensor>& field) const {
}
//------------------------------------------------------------------------------


}
