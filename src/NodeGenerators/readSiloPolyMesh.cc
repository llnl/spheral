//------------------------------------------------------------------------------
// Helper method for the SiloPolyMeshGenerator which reads a polyhedral mesh 
// from a silo file and returns the geometry.
//------------------------------------------------------------------------------
#include "Utilities/DBC.hh"
#include "Geometry/Dimension.hh"
#include "Distributed/Communicator.hh"

#include "silo.h"

#include <vector>
#include <algorithm>
using std::vector;
using std::string;
using std::pair;
using std::make_pair;

namespace Spheral {

void
readSiloPolyMesh(const std::string& fileName,
                 const std::string& meshName,
                 std::vector<Dim<3>::Vector>& positions,
                 std::vector<double>& volumes,
                 std::vector<Dim<3>::SymTensor>& H) {

  typedef Dim<3>::Vector Vector;

  // Read the mesh from the silo file.
  VERIFY2(DBInqFile(fileName.c_str()) > 0, "Error: " << fileName << " does not appear to be a valid silo file.");
  DBfile* filePtr = DBOpen(fileName.c_str(), DB_HDF5, DB_READ);
  VERIFY2(filePtr != NULL, "Error: unable to open file " << fileName);
  DBucdmesh* meshPtr = DBGetUcdmesh(filePtr, meshName.c_str());
  DBClose(filePtr);
  VERIFY2(meshPtr != NULL, "Error: unable to read mesh " << meshName << " from file " << fileName);

  // Do some checks about the mesh.
  VERIFY2(meshPtr->ndims == 3, "Error: wrong number of dimensions in mesh -- " << meshPtr->ndims);
  VERIFY2(meshPtr->phzones != NULL, "Error: missing expected phzones in mesh -- " << meshPtr->phzones);
  const unsigned nnodes = meshPtr->nnodes;
  const unsigned nzones = meshPtr->phzones->nzones;
  const unsigned nfaces = meshPtr->phzones->nfaces;

  // Compute the coordinates of the nodes.
  vector<Vector> nodePos(nnodes);
  for (unsigned i = 0; i != nnodes; ++i) {
    nodePos[i] = Vector(((double**)(meshPtr->coords))[0][i],
                        ((double**)(meshPtr->coords))[1][i],
                        ((double**)(meshPtr->coords))[2][i]);
  }

  // Find the positions and areas of all faces.
  vector<Vector> facePos(nfaces);
  vector<double> faceArea(nfaces, 0.0);
  vector<Vector> faceNormal(nfaces);
  unsigned offset = 0;
  for (unsigned i = 0; i != nfaces; ++i) {
    const unsigned n = meshPtr->phzones->nodecnt[i];
    const Vector pos1 = nodePos[meshPtr->phzones->nodelist[offset]];
    const Vector pos2 = nodePos[meshPtr->phzones->nodelist[offset + 1]];
    const Vector pos3 = nodePos[meshPtr->phzones->nodelist[offset + 2]];
    faceNormal[i] = ((pos2 - pos1).cross(pos3 - pos1)).unitVector();
    for (unsigned j = 0; j != n; ++j) {
      const unsigned jnode = meshPtr->phzones->nodelist[offset + j];
      const unsigned knode = meshPtr->phzones->nodelist[offset + ((j + 1) % n)];
      CHECK(jnode < nnodes and knode < nnodes);
      facePos[i] += nodePos[jnode];
      faceArea[i] += ((nodePos[jnode] - pos1).cross(nodePos[knode] - pos1)).magnitude();
    }
    facePos[i] /= n;
    faceArea[i] *= 0.5;
    offset += n;
  }

  // Convert the polyhedra to positions, volumes, and H tensors.
  // For now we assume tets internally, though this may be relaxed later.
  positions.resize(nzones);
  volumes.resize(nzones);
  H.resize(nzones);
  offset = 0;
  for (unsigned i = 0; i != nzones; ++i) {
    const unsigned n = meshPtr->phzones->facecnt[i];

    // First pass: find the zone position.
    for (unsigned j = 0; j != n; ++j) {
      const unsigned jface = (meshPtr->phzones->facelist[offset + j] < 0 ?
                              ~(meshPtr->phzones->facelist[offset + j]) :
                              meshPtr->phzones->facelist[offset + j]);
      CHECK(jface < nfaces);
      positions[i] += facePos[jface];
    }
    positions[i] /= n;

    // Second pass: find the volume and H.
    for (unsigned j = 0; j != n; ++j) {
      const unsigned jface = (meshPtr->phzones->facelist[offset + j] < 0 ?
                              ~(meshPtr->phzones->facelist[offset + j]) :
                              meshPtr->phzones->facelist[offset + j]);
      CHECK(jface < nfaces);
      volumes[i] += faceArea[jface] * abs((positions[i] - facePos[jface]).dot(faceNormal[jface]));
      H[i] += (facePos[jface] - positions[i]).selfdyad();
    }
    volumes[i] /= 3.0;
    H[i] = H[i].sqrt();
    CHECK(H[i].Determinant() > 0.0);
    H[i] *= pow(volumes[i]/(4.0/3.0*M_PI*H[i].Determinant()), 1.0/3.0);
    H[i] = H[i].Inverse();
    
    offset += n;
  }
}

}
