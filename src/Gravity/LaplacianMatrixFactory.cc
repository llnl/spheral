//----------------------------------------------------------------------------
//! \author $Author: jeffjohnson $
//! \date $Date: 2002-11-05 23:24:26 -0800 (Tue, 05 Nov 2002) $
//! \version $Revision: 461 $
//
// Implementation of the LaplacianMatrixFactory class.
//----------------------------------------------------------------------------

namespace Spheral {

//----------------------------------------------------------------------------
template <typename Dimension>
double** 
LaplacianMatrixFactory<Dimension>::
SPHGradSquared(const DataBase<Dimension>& db)
{
  // Count up some particles.
  const auto numberOfParticles = db.numInternalNodes();
  
  // First, we need to allocate storage for the grad and gradSquared matrices.
  // They should be N x N, where N is the number of particles.
  double** grad = new double*[numberOfParticles];
  double** gradSquared = new double*[numberOfParticles];
  for (size_t i = 0; i < numberOfParticles)
  {
    grad[i] = new double[numberOfParticles];
    std::fill(grad[i], grad[i] + numberOfParticles, 0.0);
    gradSquared[i] = new double[numberOfParticles];
    std::fill(gradSquared[i], gradSquared[i] + numberOfParticles, 0.0);
  } // end for

  // Now we've got to build the damn grad operator.

  // Some important field lists.
  const auto position = db.globalPosition();
  const auto Hfield = db.globalHfield();
  const auto& connectivityMap = db.connectivityMap();
  const auto numNodeLists = position.numFields();

  // Map each (NodeList, node) pair to the corresponding matrix index.
  std::vector<size_t> nodeListOffsets(numNodeLists + 1u, 0u);
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    nodeListOffsets[nodeListi + 1u] = (nodeListOffsets[nodeListi] +
                                      position[nodeListi]->numInternalElements());
  }
  CHECK(nodeListOffsets.back() == numberOfParticles);
  
  // Loop over each particle...
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto ni = position[nodeListi]->numInternalElements();
    for (auto nodei = 0u; nodei < ni; ++nodei) {
      const auto i = nodeListOffsets[nodeListi] + nodei;
      const auto& ri = position(nodeListi, nodei);
      const auto& Hi = Hfield(nodeListi, nodei);
      const auto& connectivity = connectivityMap.connectivityForNode(nodeListi, nodei);

      // Loop over the internal neighbors of this particle.
      for (auto nodeListj = 0u; nodeListj < numNodeLists; ++nodeListj) {
        for (const auto nodej: connectivity[nodeListj]) {
          if (nodej >= position[nodeListj]->numInternalElements()) continue;
          const auto j = nodeListOffsets[nodeListj] + nodej;
          const auto& rj = position(nodeListj, nodej);
          const auto& Hj = Hfield(nodeListj, nodej);

          // Calculate rij, the displacement between particles i and j.
          const auto rij = ri - rj;

          // Calculate the normalized displacements.
          const auto etai = Hi * rij;
          const auto etaiNorm = etai.unitVector();
          const auto etaj = Hj * rij;
          const auto etajNorm = etaj.unitVector();

          grad[i][j] = kernel.grad(r, h);
        }
      }
    }
  }

  // Now we square the grad matrix to obtain gradSquared.
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    for (size_t j = 0; j < numberOfParticles; ++j)
    {
      for (size_t k = 0; k < numberOfParticles; ++k)
      {
        gradSquared[i][j] += grad[i][k] * grad[k][j];
      } // end for
    } // end for
  } // end for

  // Toss the grad matrix.
  for (size_t i = 0; i < numberOfParticles; ++i)
  {
    delete [] grad[i];
  } // end for
  delete [] grad;
  
  // Now return the gradSquared matrix.
  return gradSquared;
} // end SPHGradSquared
//----------------------------------------------------------------------------

} // end Spheral
