//---------------------------------Spheral++----------------------------------//
// SPHGravity implementation.
//
//! \author $Author: jeffjohnson $
//! \version $Revision: 2239 $
//! \date $Date: 2007-05-28 23:58:39 -0700 (Mon, 28 May 2007) $
//----------------------------------------------------------------------------//
#include "SPHGravity.hh"
#include "Geometry/Dimension.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/NonDynamicState.hh"
#include "Hydro/WeightPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Field/Field.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Boundary/PeriodicBoundary.hh"
#include "Distributed/DistributedBoundary.hh"
#include "Utilities/DBC.hh"
#include "Material/PhysicalConstants.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Spasmos/Config.h"
#include "Spasmos/PyPetsc.h"
#include "Spasmos/MatFactory.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <algorithm>

namespace Spheral {


//------------------------------------------------------------------------------
template <typename Dimension>
SPHGravity<Dimension>::
SPHGravity(const TableKernel<Dimension>& kernel,
           typename SPHGravity<Dimension>::Scalar G,
           typename SPHGravity<Dimension>::Scalar maxDeltaVelocity /* = 2.0 */,
           typename SPHGravity<Dimension>::Scalar safetyFactor /* = 0.5 */):
  mResNorm(0.0),
  mNumIters(0),
  mG(G),
  mKernel(kernel),
  mMaxVChangeFactor(maxDeltaVelocity),
  mMinViOverAi(DBL_MAX),
  mSafetyFactor(safetyFactor),
  mMinDynTimeScale(DBL_MAX),
  mVecFactory(0),
  mMatFactory(0),
  mPotential(FieldList<Dimension, Scalar>::Copy),
  mMatrix(0),
  mRHS(0),
  mIsPeriodicOrDistributedNode(),
  mSolver(),
  mExtraEnergy(0.0),
  mNodeIndices(FieldList<Dimension, int>::Copy),
  mOverlapNodes() { 

  REQUIRE((safetyFactor > 0.0) && (safetyFactor <= 1.0));


  // Make sure Spasmos is around.
  importConfig();
  Spasmos_Initialize(SPASMOS_ERRORHANDLER_EXCEPTION);
  importPyPetsc();

  // Initialize matrix and vector factories.
  mVecFactory = VecFactory_New();
  mMatFactory = MatFactory_New();

  // Initialize a PCG solver using the normal equations.
  KSPCreate(PETSC_COMM_WORLD, &mSolver);
  KSPSetType(mSolver, KSPCGNE);

  // Use a Jacobi preconditioner.
  PC precond;
  KSPGetPC(mSolver, &precond);
  PCSetType(precond, PCJACOBI);
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
SPHGravity<Dimension>::
~SPHGravity() {
   // Get rid of the Laplacian matrix.
   Py_XDECREF(mMatrix);

   // Kill the linear solver.
   KSPDestroy(mSolver);

   // Kill the vector and matrix factories.
   Py_DECREF(mVecFactory);
   Py_DECREF(mMatFactory);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template <typename Dimension>
void 
SPHGravity<Dimension>::
mComputeMatrixStructure(const DataBase<Dimension>& dataBase,
                        const State<Dimension>& state) const
{
  // Initialize the matrix with a nonzero structure.

  // Refresh our global indexing scheme for all the nodes in the database.
  mNodeIndices = globalNodeIDs(dataBase);
  int n = mNodeIndices.numInternalNodes();

  // Make sure that the ghost nodes have the right global indexing.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mNodeIndices);
  }

  // This is necessary to ease the pain of the parallel boundary condition.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->finalizeGhostBoundary();
  }

  // Additionally, we mark ghost nodes associated with parallel and periodic 
  // boundary conditions as being special.
  mIsPeriodicOrDistributedNode.clear();
  const vector<FluidNodeList<Dimension>*>& nodeLists = dataBase.fluidNodeListPtrs();
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    if ((dynamic_cast<DistributedBoundary<Dimension>* const>(*boundaryItr) != 0) or 
        (dynamic_cast<PeriodicBoundary<Dimension>* const>(*boundaryItr) != 0))
    {
       // Get the ghost nodes for this Boundary for each node list.
       for (size_t i = 0; i < nodeLists.size(); ++i)
       {
          const vector<int>& ghostNodes = (*boundaryItr)->ghostNodes(*(nodeLists[i]));
          for (size_t j = 0; j < ghostNodes.size(); ++j)
          {
             // Mark this node as special.
             pair<int, int> nodeEntry(i, ghostNodes[j]);
             mIsPeriodicOrDistributedNode[nodeEntry] = 1;
          }
       }
    }
  }

  // Now create an "adjacency table" that contains the nonzero columns in each 
  // row for the Laplacian matrix.
  vector<vector<Py_ssize_t> > adjTable(n);
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  int numNodeLists = static_cast<int>(nodeLists.size());
  for (size_t iNodeList = 0; iNodeList < numNodeLists; ++iNodeList)
  {
    // Traverse the internal nodes in this node list.
    for (size_t iNode = 0; iNode < nodeLists[iNodeList]->numInternalNodes(); ++iNode)
    {
      // Keep track of the rows and columns in the non-zero block we create for 
      // the ith node, and automatically include the diagonal (ii) element.
      Py_ssize_t iNodeIndex = (*mNodeIndices[iNodeList])[iNode];
      Py_ssize_t row = iNodeIndex;
      adjTable[iNode].push_back(iNodeIndex);
       
      // Get the neighbors for this node.
      const vector<vector<int> >& jNodes = 
        connectivityMap.connectivityForNode(nodeLists[iNodeList], iNode);

      // Iterate over the nodelists containing the ith node's neighbors.
      for (int jNodeList = 0; jNodeList < numNodeLists; ++jNodeList)
      {
        // Get the neighboring nodes within this node list.
        const vector<int>& neighbors = jNodes[jNodeList];

        // Get the number of internal nodes within the jth nodelist.
        int numInternalNodes = nodeLists[jNodeList]->numInternalNodes();

        // Traverse these neighboring nodes.
        for (vector<int>::const_iterator jNodeIter = neighbors.begin();
             jNodeIter != neighbors.end(); ++jNodeIter)
        {
          int jNode = *jNodeIter;

          // Scrutinize non-internal nodes.  If a non-internal node is a 
          // boundary node related to a periodic or a distributed boundary, 
          // we add its contribution to the matrix.  Otherwise we discard it.
          if (jNode >= numInternalNodes)
          {
            pair<int, int> nodeEntry(jNodeList, jNode);
            if (!mIsPeriodicOrDistributedNode[nodeEntry])
              continue;
          }

          // Add the jth node to the list of rows and columns in our block.
          int jNodeIndex = ((*(mNodeIndices[jNodeList]))[jNode]);
#if 0
if ((jNodeIndex - (*mNodeIndices[0])[0] >= mNodeIndices.numInternalNodes()) ||
    (jNodeIndex < (*mNodeIndices[0])[0]))
   printf("%d: %d is off-proc\n", iNodeIndex, jNodeIndex);
#endif
          adjTable[iNode].push_back(static_cast<Py_ssize_t>(jNodeIndex));
        }
      }
    }
  }

  // Build the laplacian matrix from the given non-zero structure.
  Py_XDECREF(mMatrix);
  mMatrix = MatFactory_SparseMatrix_FromTable(mMatFactory, adjTable);

  // Tell the matrix not to allow insertions outside of the prescribed non-zeros.
  // FIXME: Something is wrong!  I don't think the nonzeros are set properly 
  // FIXME: in the parallel case!
//  MatSetOption(PyPetsc_MAT(mMatrix), MAT_NEW_NONZERO_LOCATION_ERR);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template <typename Dimension>
void 
SPHGravity<Dimension>::
mUpdateLaplacianMatrix(const DataBase<Dimension>& dataBase,
                       const State<Dimension>& state) const
{
  REQUIRE(mMatrix != 0);

  // Get the FieldLists storing the data we will use in computing our matrix.
  const FieldList<Dimension, Scalar> weight = state.scalarFields(HydroFieldNames::weight);
  const FieldList<Dimension, Vector> position = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dimension, SymTensor> Htensor = state.symTensorFields(HydroFieldNames::H);
  const vector<FluidNodeList<Dimension>*>& nodeLists = dataBase.fluidNodeListPtrs();
  int numNodeLists = static_cast<int>(nodeLists.size());
  CHECK(weight.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(Htensor.size() == numNodeLists);

  // Firstly, zero all the entries in the matrix.
  MatZeroEntries(PyPetsc_MAT(mMatrix));

  // Clear any overlapping node entries.
  mOverlapNodes.clear();

  // Traverse the list of nodes in our database and add their contributions 
  // to the matrix.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  for (size_t iNodeList = 0; iNodeList < numNodeLists; ++iNodeList)
  {
    const Field<Dimension, Scalar>& weighti = *weight[iNodeList];
    const Field<Dimension, Vector>& positioni = *position[iNodeList];
    const Field<Dimension, SymTensor>& Htensori = *Htensor[iNodeList];

    // Traverse the internal nodes in this node list.
    for (size_t iNode = 0; iNode < nodeLists[iNodeList]->numInternalNodes(); ++iNode)
    {
      // If the ith node is one that overlaps another node, we skip its 
      // row in the matrix.
      PetscInt row = static_cast<PetscInt>((*(mNodeIndices[iNodeList]))[iNode]);
      if (mOverlapNodes.find(row) != mOverlapNodes.end())
//{
//   cout << "Skipping row " << row << endl;
         continue;
//}

      // Get data for the ith node.
      Scalar wi = weighti[iNode];
      Vector xi = positioni[iNode];
      SymTensor Hi = Htensori[iNode];
      Scalar Hdeti = Hi.Determinant();

      // Keep track of the rows and columns in the non-zero block we create for 
      // the ith node, and automatically include the diagonal (ii) element.
      vector<PetscInt> columns;
      vector<double> values;
      
      // The diagonal entry should be the first one in the row, regardless of how 
      // things are actually ordered.
      columns.push_back(row);
      values.push_back(0.0);
       
      // Get the neighbors for this node.
      const vector<vector<int> >& jNodes = 
        connectivityMap.connectivityForNode(nodeLists[iNodeList], iNode);

      // Iterate over the nodelists containing the ith node's neighbors.
      for (int jNodeList = 0; jNodeList < numNodeLists; ++jNodeList)
      {
        // Get the neighboring nodes within this node list.
        const vector<int>& neighbors = jNodes[jNodeList];
        const Field<Dimension, Scalar>& weightj = *weight[jNodeList];
        const Field<Dimension, Vector>& positionj = *position[jNodeList];
        const Field<Dimension, SymTensor>& Htensorj = *Htensor[jNodeList];

        // Traverse these neighboring nodes.
        for (vector<int>::const_iterator jNodeIter = neighbors.begin();
             jNodeIter != neighbors.end(); ++jNodeIter)
        {
          int jNode = *jNodeIter;

          // Scrutinize non-internal nodes.  If a non-internal node is a 
          // boundary node related to a periodic or a distributed boundary, 
          // we add its contribution to the matrix.  Otherwise we discard it.
          if (jNode >= nodeLists[jNodeList]->numInternalNodes())
          {
            pair<int, int> nodeEntry(jNodeList, jNode);
            if (!mIsPeriodicOrDistributedNode[nodeEntry])
              continue;
          }

          // Get data for the node.
          Vector xj = positionj[jNode];
          Scalar wj = weightj[jNode];
          SymTensor Hj = Htensorj[jNode];
          Scalar Hdetj = Hj.Determinant();

          // Find the global index for this neighbor node.
          PetscInt jNodeIndex = 
            static_cast<PetscInt>((*(mNodeIndices[jNodeList]))[jNode]);

          // Node displacement.
          const Vector xij = xi - xj;
          const Vector etai = Hi*xij;
          const Vector etaj = Hj*xij;
          const Scalar etaMagi = etai.magnitude();
          const Scalar etaMagj = etaj.magnitude();
          CHECK(etaMagi >= 0.0);
          CHECK(etaMagj >= 0.0);

          // Does node j overlap node i?  That is, are the two so close 
          // together that they are physically indistinguishable?  If so,
          // skip this entry of the matrix.
          if (etaMagi < 0.1)
          {
             mOverlapNodes[jNodeIndex] = row;
//printf("Overlapping node %d -> %d\n", jNodeIndex, row);
             continue;
          }

          // Add the jth node to the list of columns.
          columns.push_back(jNodeIndex);

          // Symmetrized kernel weight and gradient.
          const Scalar Wi = mKernel(etaMagi, Hdeti);
          const Vector Hetai = Hi*etai.unitVector();
          const Vector gradWi = Hetai*mKernel.grad(etaMagi, Hdeti);

          const Scalar Wj = mKernel(etaMagj, Hdetj);
          const Vector Hetaj = Hj*etaj.unitVector();
          const Vector gradWj = Hetaj*mKernel.grad(etaMagj, Hdetj);

          const Vector gradWij = 0.5 * (gradWi + gradWj);

          // Now compute the off-diagonal value, softening the interaction 
          // of nodes that are close together.
          double hi = Dimension::rootnu(Hdeti);
          double Rij = xij.magnitude2() + 0.01*hi*hi;
          double offDiagValue = 2.0 * wj * gradWij.dot(xij)/Rij;
#if 0
if ((jNodeIndex - (*mNodeIndices[0])[0] >= mNodeIndices.numInternalNodes()) ||
    (jNodeIndex < (*mNodeIndices[0])[0]))
   printf("%d: Off-proc contrib from %d is %g\n", row, jNodeIndex, offDiagValue);
//cout << "[" << row << "," <<  jNodeIndex << "]: " << offDiagValue << " (|xij| = " << xij.magnitude2() << ")" << endl;
#endif
          values.push_back(offDiagValue);

          // The diagonal value for this row is the negative sum of the off-diagonal
          // values, so we simply subtract.
          values[0] -= offDiagValue;
        }
      }
      assert(values.size() == columns.size());

      // Now insert the row into the matrix. 
      if (MatSetValues(PyPetsc_MAT(mMatrix), 1, &row, 
                       static_cast<PetscInt>(columns.size()), &(columns[0]),
                       &(values[0]), INSERT_VALUES) != 0)
      {
         PyErr_Print();
         VERIFY(false);
      }
    }
  }

  // Before we leave, we must assemble the matrix.
  MatAssemblyBegin(PyPetsc_MAT(mMatrix), MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(PyPetsc_MAT(mMatrix), MAT_FINAL_ASSEMBLY);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template <typename Dimension>
PyObject* 
SPHGravity<Dimension>::
mCreateRHS(const State<Dimension>& state) const {

  REQUIRE(mMatrix != 0);
//  REQUIRE(PetscMatrix_Check(mMatrix));

  // Here's a brand new vector that fits with our matrix!
  Py_ssize_t n = static_cast<Py_ssize_t>(mPotential.numInternalNodes());
  PyObject* RHS = VecFactory_Vector(mVecFactory, n);
  CHECK(RHS != 0);

  // We go over the interior nodes and stick the mass density in there.
  const FieldList<Dimension, Scalar> rho = state.scalarFields(HydroFieldNames::massDensity);
  vector<double> localVec(rho.numInternalNodes());
  size_t offset = 0;
  for (int iNodeList = 0; iNodeList < rho.numFields(); ++iNodeList) {
    Field<Dimension, Scalar>& rhoi = *(rho[iNodeList]);
    int numInternalNodes = rhoi.nodeList().numInternalNodes();
    for (int iNode = 0; iNode < numInternalNodes; ++iNode) {
      localVec[offset] = -4.0*M_PI*mG*rhoi[iNode];
      ++offset;
    }
  }
  assert(offset == rho.numInternalNodes());

  // Add contributions from boundary nodes to the RHS.
  const vector<Boundary<Dimension>*>& BCs = this->boundaryConditions();
  for (size_t i = 0; i < BCs.size(); ++i)
  {
     // Ignore periodic boundary nodes, since we've already accounted for 
     // them in the matrix.
    if ((dynamic_cast<PeriodicBoundary<Dimension>* const>(BCs[i]) == 0) &&
        (dynamic_cast<DistributedBoundary<Dimension>* const>(BCs[i]) == 0))
     {
        // FIXME: Fun stuff goes here.
     }
  }

  // Splat the values back into the RHS vector.
  PyVec_CopyFromArray(RHS, 0, static_cast<Py_ssize_t>(localVec.size()),
                      &(localVec[0]));

  // Now go over the set of overlapping nodes and mash their mass densities
  // together.
  if (!mOverlapNodes.empty())
  {
    vector<PetscInt> overlappingRows, degenerateRows;
    for (map<int, int>::const_iterator iter = mOverlapNodes.begin();
         iter != mOverlapNodes.end(); ++iter)
    {
      degenerateRows.push_back(static_cast<PetscInt>(iter->first));
      overlappingRows.push_back(static_cast<PetscInt>(iter->second));
    }
    vector<double> degenerateDensities(mOverlapNodes.size());
    VecGetValues(PyPetsc_VEC(RHS), static_cast<PetscInt>(degenerateRows.size()),
          &(degenerateRows[0]), &(degenerateDensities[0]));
    VecSetValues(PyPetsc_VEC(RHS), static_cast<PetscInt>(overlappingRows.size()),
          &(overlappingRows[0]), &(degenerateDensities[0]), ADD_VALUES);
    VecAssemblyBegin(PyPetsc_VEC(RHS));
    VecAssemblyEnd(PyPetsc_VEC(RHS));

    // Make sure to zero out the entries on the degenerate rows, too. 
    fill(degenerateDensities.begin(), degenerateDensities.end(), 0.0);
    VecSetValues(PyPetsc_VEC(RHS), static_cast<PetscInt>(degenerateRows.size()),
          &(degenerateRows[0]), &(degenerateDensities[0]), INSERT_VALUES);
    VecAssemblyBegin(PyPetsc_VEC(RHS));
    VecAssemblyEnd(PyPetsc_VEC(RHS));
  }

  // Outta here.
  return RHS;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template <typename Dimension>
void
SPHGravity<Dimension>::
mComputeGravitationalPotential(const DataBase<Dimension>& dataBase,
                               const State<Dimension>& state) const {
  // First, make sure that the matrix is up to date. 
  REQUIRE(mMatrix != 0);

  // Create a vector representing the RHS of the Poisson equation.
  Py_XDECREF(mRHS);
  mRHS = mCreateRHS(state);
  CHECK(mRHS != 0);

  // Now set up the linear system and solve it.
  int n = mPotential.numInternalNodes();
  PyObject* sol = VecFactory_Vector(mVecFactory, n);
  KSPSetOperators(mSolver, PyPetsc_MAT(mMatrix), 
                  PyPetsc_MAT(mMatrix), SAME_PRECONDITIONER);
  KSPSolve(mSolver, PyPetsc_VEC(mRHS), PyPetsc_VEC(sol));

  // Get data about the linear solution.
  KSPGetResidualNorm(mSolver, &mResNorm);
  KSPGetIterationNumber(mSolver, &mNumIters);

  // Did the solution converge?  If not, we have to let the user know What 
  // Happened.  The easiest way to do this is to let Spasmos decipher the 
  // PETSc error.
  if (!PyPetsc_SolveConverged((void*)mSolver, "SPH Gravity Poisson equation"))
  {
    PyErr_Print();
    VERIFY(false);
  }

  // Now go over the set of overlapping nodes and copy the solutions back 
  // to the associated degenerate nodes.
  // together.
  if (!mOverlapNodes.empty())
  {
    vector<PetscInt> overlappingRows, degenerateRows;
    for (map<int, int>::const_iterator iter = mOverlapNodes.begin();
         iter != mOverlapNodes.end(); ++iter)
    {
      degenerateRows.push_back(static_cast<PetscInt>(iter->first));
      overlappingRows.push_back(static_cast<PetscInt>(iter->second));
//printf("Copying solution from %d to %d\n", iter->second, iter->first);
    }
    vector<double> degenerateDensities(mOverlapNodes.size());
    VecGetValues(PyPetsc_VEC(sol), static_cast<PetscInt>(overlappingRows.size()),
                 &(overlappingRows[0]), &(degenerateDensities[0]));
    VecSetValues(PyPetsc_VEC(sol), static_cast<PetscInt>(degenerateRows.size()),
                 &(degenerateRows[0]), &(degenerateDensities[0]), INSERT_VALUES);
    VecAssemblyBegin(PyPetsc_VEC(sol));
    VecAssemblyEnd(PyPetsc_VEC(sol));
  }

  // Transfer the data in the solution to our gravitational potential field 
  // and find its maximum value.
//  CHECK(mPotential.numInternalNodes() == PetscVector_LocalSize(sol));
  vector<double> localVec(mPotential.numInternalNodes());
  double maxValue = -FLT_MAX;
  PyVec_CopyToArray(sol, 0, static_cast<Py_ssize_t>(localVec.size()),
                    &(localVec[0]));
  for (int iNodeList = 0; iNodeList < mPotential.numFields(); ++iNodeList) {
    Field<Dimension, Scalar>& potential = *(mPotential[iNodeList]);
    int numInternalNodes = potential.nodeList().numInternalNodes();
    for (int iNode = 0; iNode < numInternalNodes; ++iNode) {
      potential[iNode] = localVec[iNode];
      if (potential[iNode] > maxValue)
         maxValue = potential[iNode];
    }
  }

  // Now subtract the maximum value from the potential to give us a 
  // well-defined potential with a maximum value of zero.
  for (int iNodeList = 0; iNodeList < mPotential.numFields(); ++iNodeList) {
    Field<Dimension, Scalar>& potential = *(mPotential[iNodeList]);
    int numInternalNodes = potential.nodeList().numInternalNodes();
    for (int iNode = 0; iNode < numInternalNodes; ++iNode) {
      potential[iNode] -= maxValue;
    }
  }

  // Clean up.
  Py_DECREF(sol);
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template <typename Dimension>
void 
SPHGravity<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const
{

  // Access to pertinent fields in the database.
  const FieldList<Dimension, Scalar> m = state.scalarFields(HydroFieldNames::mass);
  const FieldList<Dimension, Scalar> rho = state.scalarFields(HydroFieldNames::massDensity);
  const FieldList<Dimension, Vector> x = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dimension, Vector> v = state.vectorFields(HydroFieldNames::velocity);
  const FieldList<Dimension, SymTensor> H = state.symTensorFields(HydroFieldNames::H);

  // Get the accelerations we'll be modifying.
  FieldList<Dimension, Vector> DvDt = derivs.vectorFields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity);

  // Zero out the total gravitational potential energy.
  mExtraEnergy = 0.0;

  // Update the Laplacian matrix.
  mUpdateLaplacianMatrix(dataBase, state);

  // Compute the gravitational potential.
  mComputeGravitationalPotential(dataBase, state);

  // Make sure that the potential is correct at the boundaries.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mPotential);
  }
  // This is necessary to ease the pain of the parallel boundary condition.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->finalizeGhostBoundary();
  }


  // Now compute the gravitational acceleration, which is the negative gradient of the 
  // potential, and add it to the nodal accelerations.
  mMinViOverAi = DBL_MAX;
  mMinDynTimeScale = DBL_MAX;
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  int numNodeLists = DvDt.numFields();
  for (int iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
    const NodeList<Dimension>& nodeList = DvDt[iNodeList]->nodeList(); 
    for (int iNode = 0; iNode < nodeList.numInternalNodes(); ++iNode)
    {
      // Get a reference to the acceleration vector and the gravitational potential.
      Vector& ai = (*DvDt[iNodeList])[iNode];
      const SymTensor& Hi = (*H[iNodeList])[iNode];
      Scalar detHi = Hi.Determinant();
      Scalar phii = (*mPotential[iNodeList])[iNode];
      Scalar rhoi = (*rho[iNodeList])[iNode];
    
      // Loop over the neighboring nodes.
      const vector< vector<int> >& fullConnectivity = 
         connectivityMap.connectivityForNode(&nodeList, iNode);
      // Iterate over the NodeLists.
      for (int jNodeList = 0; jNodeList != numNodeLists; ++jNodeList) {
        const vector<int>& connectivity = fullConnectivity[jNodeList];
        if (connectivity.size() > 0) {
          // Loop over the neighbors.
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            int jNode = *jItr;
            if (iNode != jNode) { // No self-interaction. 
              // Contribute to the gradient of the gravitational potential.
              Scalar phij = (*mPotential[jNodeList])[jNode];
              Scalar mj = (*m[jNodeList])[jNode];
              Vector xij = (*x[iNodeList])[iNode] - 
                           (*x[jNodeList])[jNode];
              CHECK(xij.magnitude2() != 0.0);

              const SymTensor& Hj = (*H[jNodeList])[jNode];
              Scalar detHj = Hj.Determinant();

              Vector etai = Hi.dot(xij);
              Scalar dWi = mKernel.gradValue(etai.magnitude(), detHi);
              Vector etaiHat = etai.unitVector();
              Vector gradWi = dWi * Hi.dot(etaiHat);

              Vector etaj = Hj.dot(xij);
              Scalar dWj = mKernel.gradValue(etaj.magnitude(), detHj);
              Vector etajHat = etaj.unitVector();
              Vector gradWj = dWj * Hj.dot(etajHat);

              // Use the form of the SPH gradient that shoves the density 
              // into the gradient.
              Vector gradWij = 0.5 * (gradWi + gradWj);
              ai -= mj * (phij-phii) * gradWij / rhoi;
            }
          }
        }
      } 
      // Add this node's contribution to the potential energy.
      mExtraEnergy += (*m[iNodeList])[iNode] * (*mPotential[iNodeList])[iNode];

      // Capture the minimum velocity/acceleration ratio.
      if (ai.magnitude() > 0.0)
      {
         Scalar vtoaRatio = (*v[iNodeList])[iNode].magnitude()/ai.magnitude();
         if (mMinViOverAi > vtoaRatio)
            mMinViOverAi = vtoaRatio;
      }

      // Capture the minimum dynamical timescale.
      if (rhoi > 0.0)
      {
         Scalar minDynTimeScale = 1.0/sqrt(mG*rhoi);
         if (mMinDynTimeScale > minDynTimeScale)
            mMinDynTimeScale = minDynTimeScale;
      }
    } 
  } // end for

  // That's it.
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHGravity<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::ScalarPolicyPointerType ScalarPolicyPointer;
  typedef typename State<Dimension>::VectorPolicyPointerType VectorPolicyPointer;
  typedef typename State<Dimension>::TensorPolicyPointerType TensorPolicyPointer;
  typedef typename State<Dimension>::SymTensorPolicyPointerType SymTensorPolicyPointer;

  // Allocate space for the gravitational potential FieldList if necessary.
  if (mPotential.numFields() == 0)
  {
    mPotential = dataBase.newGlobalFieldList(Scalar());
  } // end if

  // Traverse the set of fluid node lists in the database.
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {

    // Stuff used to figure out the gravitational potential. 
    ScalarPolicyPointer massPolicy(new NonDynamicState<Dimension, Field<Dimension, Scalar> >());
    ScalarPolicyPointer weightPolicy(new WeightPolicy<Dimension>());
    ScalarPolicyPointer rhoPolicy(new NonDynamicState<Dimension, Field<Dimension, Scalar> >());
    VectorPolicyPointer xPolicy(new NonDynamicState<Dimension, Field<Dimension, Vector> >());
    SymTensorPolicyPointer HPolicy(new NonDynamicState<Dimension, Field<Dimension, SymTensor> >());

    state.registerField((*itr)->mass(), massPolicy);
    state.registerField((*itr)->weight(), weightPolicy);
    state.registerField((*itr)->massDensity(), rhoPolicy);
    state.registerField((*itr)->positions(), xPolicy);
    state.registerField((*itr)->Hfield(), HPolicy);

    // Velocity update.
    VectorPolicyPointer velocityPolicy(new IncrementState<Dimension, Field<Dimension, Vector> >());
    state.registerField((*itr)->velocity(), velocityPolicy);
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHGravity<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  int i = 0;
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++i) {
    derivs.registerField((*itr)->DvelocityDt());
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHGravity<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  REQUIRE(this->valid());

#if 0
  // Apply boundary conditions to the pertinent fields.
//  FieldList<Dimension, Vector> position = state.vectorFields(HydroFieldNames::position);
  FieldList<Dimension, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
  FieldList<Dimension, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
  FieldList<Dimension, Scalar> massDensity = state.scalarFields(HydroFieldNames::massDensity);
  FieldList<Dimension, Scalar> weight = state.scalarFields(HydroFieldNames::weight);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(H);
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(weight);
  }
#endif
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template<typename Dimension>
void
SPHGravity<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  REQUIRE(this->valid());

#if 0
  // Enforce boundary conditions on the fluid state Fields.
  FieldList<Dimension, Vector> position = state.vectorFields(HydroFieldNames::position);
  FieldList<Dimension, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
  FieldList<Dimension, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
  FieldList<Dimension, Scalar> massDensity = state.scalarFields(HydroFieldNames::massDensity);
  FieldList<Dimension, Scalar> weight = state.scalarFields(HydroFieldNames::weight);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(position);
    (*boundaryItr)->enforceFieldListBoundary(H);
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(weight);
  }
#endif
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
bool
SPHGravity<Dimension>::
initialize(const Scalar& time, 
           const Scalar& dt,
           const DataBase<Dimension>& db, 
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs)
{
  mComputeMatrixStructure(db, state);
  return false;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
typename SPHGravity<Dimension>::TimeStepType
SPHGravity<Dimension>::
dt(const DataBase<Dimension>& dataBase, 
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const
{
  // The maximum change in our velocity dictates the maximum allowable time step.
  // If ai is the acceleration experienced by an SPH node, and vi is that node's 
  // velocity, the factor by which its velocity changes over a time interval dti is 
  // abs(ai)/abs(vi) * dti.  If the maximum allowed velocity change factor is F, then 
  // the maximum allowable time step is dt = F * min(abs(vi)/abs(ai)). 

  // If this is the first time step, we won't have computed the minimum (v/a) yet, so we'll do 
  // it here.  We can estimate the maximum allowed time step from the gravitational acceleration, 
  // assuming that any hydro accelerations are accounted for in the SPHGravity package's 
  // bookkeeping.  This allows us to run gravity simulations on "dust-like" (P = 0) mass 
  // distributions.
  if (mMinViOverAi < 0) {
    FieldList<Dimension, Vector> velocity = state.vectorFields(HydroFieldNames::velocity);
    FieldList<Dimension, Vector> accel = state.vectorFields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity);
    mMinViOverAi = DBL_MAX;
    int numNodeLists = velocity.numFields();
    for (int iNodeList = 0; iNodeList < numNodeLists; ++iNodeList) {
      int numInternalNodes = velocity[iNodeList]->nodeList().numInternalNodes();
      const Field<Dimension, Vector>& v = *(velocity[iNodeList]);
      const Field<Dimension, Vector>& a = *(accel[iNodeList]);
      for (int iNode = 0; iNode < numInternalNodes; ++iNode) {
        Scalar viOverAi = v[iNode].magnitude() / a[iNode].magnitude();
        if (viOverAi < mMinViOverAi)
          mMinViOverAi = viOverAi;
      }
    }
  }

  const double deltatV = mMaxVChangeFactor * mMinViOverAi;
  if (deltatV < mMinDynTimeScale)
  {
    stringstream reasonStream;
    reasonStream << "Maximum velocity change factor" << endl
                 << "min(v/a): " << mMinViOverAi << endl
                 << "max V change factor: " << mMaxVChangeFactor << endl
                 << "dt = safetyFactor * max V change factor * min(v/a): " << mSafetyFactor * deltatV << ends;
    return TimeStepType(mSafetyFactor * deltatV, reasonStream.str());
  }
  else
  {
    stringstream reasonStream;
    reasonStream << "Minimum gravitational dynamic time scale" << endl
                 << "min(1/sqrt(G*rho)): " << mMinDynTimeScale << endl
                 << "dt = safetyFactor * min(1/sqrt(G*rho)): " << mSafetyFactor * mMinDynTimeScale << ends;
    return TimeStepType(mSafetyFactor * mMinDynTimeScale, reasonStream.str());
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
typename SPHGravity<Dimension>::Scalar 
SPHGravity<Dimension>::
extraEnergy() const
{
  return mExtraEnergy;
}
//------------------------------------------------------------------------------


//------------------------------------------------------------------------------
template <typename Dimension>
const FieldList<Dimension, typename SPHGravity<Dimension>::Scalar>&
SPHGravity<Dimension>::
potential() const
{
  return mPotential;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
template <typename Dimension>
bool 
SPHGravity<Dimension>::
valid() const
{
  // Make sure that we have only parallel and periodic boundary 
  // conditions for now.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    if ((dynamic_cast<DistributedBoundary<Dimension>* const>(*boundaryItr) == 0) and 
        (dynamic_cast<PeriodicBoundary<Dimension>* const>(*boundaryItr) == 0))
    {
       return false;
    }
  }
  return true;
}

//------------------------------------------------------------------------------
} // end namespace Spheral

