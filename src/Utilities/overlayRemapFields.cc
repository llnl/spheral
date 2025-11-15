//------------------------------------------------------------------------------
// Use geometric clipping to remap a set of conserved fields.
// Currently only works single NodeList -> single NodeList, no boundaries.
//------------------------------------------------------------------------------
#include "overlayRemapFields.hh"
#include "Utilities/clipFacetedVolume.hh"
#include "DataBase/DataBase.hh"
#include "VoronoiCells/computeVoronoiVolume.hh"
#include "Geometry/GeomPlane.hh"
#include "Utilities/DBC.hh"
#include "Distributed/Process.hh"
#include "Distributed/Communicator.hh"

#include <map>
#include <algorithm>
using std::vector;
using std::tuple;
using std::list;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

template<typename Dimension>
void
overlayRemapFields(const vector<Boundary<Dimension>*>& boundaries,
                   const vector<Field<Dimension, typename Dimension::Scalar>*>& scalarDonorFields,
                   const vector<Field<Dimension, typename Dimension::Vector>*>& vectorDonorFields,
                   const vector<Field<Dimension, typename Dimension::Tensor>*>& tensorDonorFields,
                   const vector<Field<Dimension, typename Dimension::SymTensor>*>& symTensorDonorFields,
                   vector<Field<Dimension, typename Dimension::Scalar>*>& scalarAcceptorFields,
                   vector<Field<Dimension, typename Dimension::Vector>*>& vectorAcceptorFields,
                   vector<Field<Dimension, typename Dimension::Tensor>*>& tensorAcceptorFields,
                   vector<Field<Dimension, typename Dimension::SymTensor>*>& symTensorAcceptorFields) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;
  typedef typename FacetedVolume::Facet Facet;
  typedef GeomPlane<Dimension> Plane;


  // Preconditions.
  VERIFY2(scalarDonorFields.size() == scalarAcceptorFields.size(), "overlayRemapFields ERROR: number of acceptor scalar fields does not match number of donors.");
  VERIFY2(vectorDonorFields.size() == vectorAcceptorFields.size(), "overlayRemapFields ERROR: number of acceptor vector fields does not match number of donors.");
  VERIFY2(tensorDonorFields.size() == tensorAcceptorFields.size(), "overlayRemapFields ERROR: number of acceptor tensor fields does not match number of donors.");
  VERIFY2(symTensorDonorFields.size() == symTensorAcceptorFields.size(), "overlayRemapFields ERROR: number of acceptor symTensor fields does not match number of donors.");
  const unsigned nScalarFields = scalarDonorFields.size();
  const unsigned nVectorFields = vectorDonorFields.size();
  const unsigned nTensorFields = tensorDonorFields.size();
  const unsigned nSymTensorFields = symTensorDonorFields.size();

  // Find the donor and acceptor NodeLists.
  const NodeList<Dimension> *donorNodeListPtr = NULL, *acceptorNodeListPtr = NULL;
  for (const Field<Dimension, Scalar>* field: scalarDonorFields) {
    VERIFY2(donorNodeListPtr == NULL or donorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all donor fields must be on same NodeList.");
    donorNodeListPtr = field->nodeListPtr();
  }
  for (const Field<Dimension, Vector>* field: vectorDonorFields) {
    VERIFY2(donorNodeListPtr == NULL or donorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all donor fields must be on same NodeList.");
    donorNodeListPtr = field->nodeListPtr();
  }
  for (const Field<Dimension, Tensor>* field: tensorDonorFields) {
    VERIFY2(donorNodeListPtr == NULL or donorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all donor fields must be on same NodeList.");
    donorNodeListPtr = field->nodeListPtr();
  }
  for (const Field<Dimension, SymTensor>* field: symTensorDonorFields) {
    VERIFY2(donorNodeListPtr == NULL or donorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all donor fields must be on same NodeList.");
    donorNodeListPtr = field->nodeListPtr();
  }
  for (Field<Dimension, Scalar>* field: scalarAcceptorFields) {
    VERIFY2(acceptorNodeListPtr == NULL or acceptorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all acceptor fields must be on same NodeList.");
    acceptorNodeListPtr = field->nodeListPtr();
    *field = 0.0;
  }
  for (Field<Dimension, Vector>* field: vectorAcceptorFields) {
    VERIFY2(acceptorNodeListPtr == NULL or acceptorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all acceptor fields must be on same NodeList.");
    acceptorNodeListPtr = field->nodeListPtr();
    *field = Vector::zero();
  }
  for (Field<Dimension, Tensor>* field: tensorAcceptorFields) {
    VERIFY2(acceptorNodeListPtr == NULL or acceptorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all acceptor fields must be on same NodeList.");
    acceptorNodeListPtr = field->nodeListPtr();
    *field = Tensor::zero();
  }
  for (Field<Dimension, SymTensor>* field: symTensorAcceptorFields) {
    VERIFY2(acceptorNodeListPtr == NULL or acceptorNodeListPtr == field->nodeListPtr(), "overlayRemapFields ERROR: all acceptor fields must be on same NodeList.");
    acceptorNodeListPtr = field->nodeListPtr();
    *field = SymTensor::zero();
  }
  const unsigned nD = donorNodeListPtr->numInternalNodes();
  Neighbor<Dimension>& neighborD = donorNodeListPtr->neighbor();
  Neighbor<Dimension>& neighborA = acceptorNodeListPtr->neighbor();

  // Build the donor volumes.
  Field<Dimension, FacetedVolume> localDonorCells("donor cells", *donorNodeListPtr);
  {
    DataBase<Dimension> db;
    db.appendNodeList(*dynamic_cast<FluidNodeList<Dimension>*>(const_cast<NodeList<Dimension>*>(donorNodeListPtr)));
    neighborD.updateNodes();
    for (Boundary<Dimension>* boundPtr: boundaries) boundPtr->setAllGhostNodes(db);
    for (Boundary<Dimension>* boundPtr: boundaries) boundPtr->finalizeGhostBoundary();
    neighborD.updateNodes();
    db.updateConnectivityMap(false, false, false);
    const auto& cm = db.connectivityMap();
    const auto position = db.fluidPosition();
    const auto H = db.fluidHfield();
    const auto damage = db.solidDamage();
    const auto weight = db.newFluidFieldList(1.0, "weight");
    auto etaVoidPoints = db.newFluidFieldList(vector<Vector>(), "eta void points");
    auto surfacePoint = db.newFluidFieldList(0, "surface point");
    auto vol = db.newFluidFieldList(0.0, "volume");
    auto deltaMedian = db.newFluidFieldList(Vector::zero(), "displacement");
    FieldList<Dimension, FacetedVolume> cells_fl(FieldStorageType::ReferenceFields);
    cells_fl.appendField(localDonorCells);
    auto cellFaceFlags_fl = db.newFluidFieldList(vector<CellFaceFlag>(), "face flags");
    computeVoronoiVolume(position, H, cm, damage, vector<FacetedVolume>(), vector<vector<FacetedVolume>>(), boundaries, weight, 
                         surfacePoint, vol, deltaMedian, etaVoidPoints, cells_fl, cellFaceFlags_fl);
    const_cast<NodeList<Dimension>*>(donorNodeListPtr)->numGhostNodes(0);
    neighborD.updateNodes();
    // cerr << "Donor volume range: " << vol.min() << " " << vol.max() << endl;
  }
    
  // Build the acceptor volumes.
  Field<Dimension, FacetedVolume> localAcceptorCells("acceptor cells", *acceptorNodeListPtr);
  {
    neighborA.updateNodes();
    DataBase<Dimension> db;
    db.appendNodeList(*dynamic_cast<FluidNodeList<Dimension>*>(const_cast<NodeList<Dimension>*>(acceptorNodeListPtr)));
    for (Boundary<Dimension>* boundPtr: boundaries) boundPtr->setAllGhostNodes(db);
    for (Boundary<Dimension>* boundPtr: boundaries) boundPtr->finalizeGhostBoundary();
    neighborA.updateNodes();
    db.updateConnectivityMap(false, false, false);
    const auto& cm = db.connectivityMap();
    const auto position = db.fluidPosition();
    const auto H = db.fluidHfield();
    const auto damage = db.solidDamage();
    const auto weight = db.newFluidFieldList(1.0, "weight");
    auto etaVoidPoints = db.newFluidFieldList(vector<Vector>(), "eta void points");
    auto surfacePoint = db.newFluidFieldList(0, "surface point");
    auto vol = db.newFluidFieldList(0.0, "volume");
    auto deltaMedian = db.newFluidFieldList(Vector::zero(), "displacement");
    FieldList<Dimension, FacetedVolume> cells_fl(FieldStorageType::ReferenceFields);
    cells_fl.appendField(localAcceptorCells);
    auto cellFaceFlags_fl = db.newFluidFieldList(vector<CellFaceFlag>(), "face flags");
    computeVoronoiVolume(position, H, cm, damage, vector<FacetedVolume>(), vector<vector<FacetedVolume>>(), boundaries, weight, 
                         surfacePoint, vol, deltaMedian, etaVoidPoints, cells_fl, cellFaceFlags_fl);
    const_cast<NodeList<Dimension>*>(acceptorNodeListPtr)->numGhostNodes(0);
    neighborA.updateNodes();
    // cerr << "Acceptor volume range: " << vol.min() << " " << vol.max() << endl;
  }

  const Field<Dimension, Vector>& posD = donorNodeListPtr->positions();
  const Field<Dimension, SymTensor>& HD = donorNodeListPtr->Hfield();
  // const Field<Dimension, Vector>& posA = acceptorNodeListPtr->positions();

#ifdef SPHERAL_ENABLE_MPI
  // Parallel info
  const int myproc = Process::getRank();
  const int nprocs = Process::getTotalNumberOfProcesses();
  //..........................................................................
  // Parallel version
  // Pack up and broadcast our donor volumes to everyone else.
  vector<char> localPackedDonorCells;
  packElement(nD, localPackedDonorCells);
  for (auto i = 0u; i != nD; ++i) packElement(localDonorCells(i), localPackedDonorCells);
  for (auto i = 0u; i != nD; ++i) packElement(posD(i), localPackedDonorCells);
  for (auto i = 0u; i != nD; ++i) packElement(HD(i), localPackedDonorCells);
  unsigned nlocalpack = localPackedDonorCells.size();

  // Prepare buffers for asynchronous sends.
  vector<unsigned> sendBufsizes;
  vector<MPI_Request> sendRequests;
  sendBufsizes.reserve(4*nprocs);
  sendRequests.reserve(4*nprocs);

  // Look for intersecting node volumes.
  vector<vector<unsigned>> intersectDonorIndices(nprocs);                 // index of donor node    [domain][donorID]
  vector<vector<vector<unsigned>>> intersectAcceptorIndices(nprocs);      // index of acceptor node [domain][donorID][acceptorID]
  vector<vector<vector<double>>>  intersectVols(nprocs);                  // volume of intersection [domain][donorID][acceptorID]
  list<vector<char>> sendBuffers;
  for (int iproc = 0; iproc != nprocs; ++iproc) {
    if (myproc == 0) cerr << "Intersecting domain " << iproc << " of " << nprocs << "..." << endl;

    // Get the other processes donor info.
    unsigned npack = nlocalpack;
    MPI_Bcast(&npack, 1, MPI_UNSIGNED, iproc, Communicator::communicator());
    vector<char> buffer(npack);
    if (iproc == myproc) buffer = localPackedDonorCells;
    MPI_Bcast(&(*buffer.begin()), npack, MPI_CHAR, iproc, Communicator::communicator());
    vector<char>::const_iterator bufItr = buffer.begin();
    unsigned donorN;
    unpackElement(donorN, bufItr, buffer.end());
    vector<FacetedVolume> donorCells(donorN);
    vector<Vector> donorPos(donorN);
    vector<SymTensor> donorH(donorN);
    for (auto i = 0u; i != donorN; ++i) unpackElement(donorCells[i], bufItr, buffer.end());
    for (auto i = 0u; i != donorN; ++i) unpackElement(donorPos[i], bufItr, buffer.end());
    for (auto i = 0u; i != donorN; ++i) unpackElement(donorH[i], bufItr, buffer.end());
    CHECK(bufItr == buffer.end());

    for (unsigned i = 0; i != donorN; ++i) {
      vector<int> Amasters, Acoarse, Arefine;
      neighborA.setMasterList(donorPos[i], donorH[i], Amasters, Acoarse);
      neighborA.setRefineNeighborList(donorPos[i], donorH[i], Acoarse, Arefine);
      for (auto jitr = Arefine.begin(); jitr < Arefine.end(); ++jitr) {
        const int j = *jitr;
        if (donorCells[i].intersect(localAcceptorCells(j))) {
          const vector<Facet>& facets = localAcceptorCells(j).facets();
          vector<Plane> planes;
          planes.reserve(facets.size());
          for (const Facet& facet: facets) planes.push_back(Plane(facet.position(), -facet.normal()));
          const Scalar Vi = clippedVolume(donorCells[i], planes);
          if (Vi > 0.0) {
            if (intersectDonorIndices[iproc].empty() or intersectDonorIndices[iproc].back() != i) {
              intersectDonorIndices[iproc].push_back(i);
              intersectAcceptorIndices[iproc].push_back(vector<unsigned>());
              intersectVols[iproc].push_back(vector<double>());
            }
            intersectAcceptorIndices[iproc].back().push_back(j);
            intersectVols[iproc].back().push_back(Vi);
          }
        }
      }
      CHECK(intersectDonorIndices[iproc].size() == intersectAcceptorIndices[iproc].size());
      CHECK(intersectDonorIndices[iproc].size() == intersectVols[iproc].size());
    }

    // Send to each process the donors we are using from them, along with the total volume intersected for each donor.
    {
      const unsigned n = intersectDonorIndices[iproc].size();
      CHECK(intersectVols[iproc].size() == n);
      sendBuffers.push_back(vector<char>());
      packElement(intersectDonorIndices[iproc], sendBuffers.back());
      for (auto i = 0u; i != n; ++i) packElement(accumulate(intersectVols[iproc][i].begin(), intersectVols[iproc][i].end(), 0.0), sendBuffers.back());
      sendBufsizes.push_back(sendBuffers.back().size());
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&sendBufsizes.back(), 1, MPI_UNSIGNED, iproc, 20, Communicator::communicator(), &sendRequests.back());
      if (sendBufsizes.back() > 0) {
        sendRequests.push_back(MPI_Request());
        MPI_Isend(&(*sendBuffers.back().begin()), sendBufsizes[iproc], MPI_CHAR, iproc, 21, Communicator::communicator(), &sendRequests.back());
      }
    }
  }

  // Get the receive information about our donors from other domains.
  // This includs the indices of our donors used by each domain (stored to donorNodeIDs), and the cumulative volume used per donor node on each domain.
  // We accumulate the total intersected volume for each donor node to the voltot field.
  vector<vector<unsigned>> donorNodeIDs(nprocs);
  Field<Dimension, Scalar> voltot("intersected volume totals", *donorNodeListPtr);
  for (int iproc = 0; iproc != nprocs; ++iproc) {
    unsigned bufsize;
    MPI_Status recvStatus;
    MPI_Recv(&bufsize, 1, MPI_UNSIGNED, iproc, 20, Communicator::communicator(), &recvStatus);
    if (bufsize > 0) {
      vector<char> buffer(bufsize);
      MPI_Status recvStatus;
      MPI_Recv(&buffer[0], bufsize, MPI_CHAR, iproc, 21, Communicator::communicator(), &recvStatus);
      vector<char>::const_iterator bufItr = buffer.begin();
      unpackElement(donorNodeIDs[iproc], bufItr, buffer.end());
      const unsigned n = donorNodeIDs[iproc].size();
      CHECK(n <= nD);
      double voli;
      for (auto i = 0u; i != n; ++i) {
        unpackElement(voli, bufItr, buffer.end());
        voltot(donorNodeIDs[iproc][i]) += voli;
      }
      CHECK(bufItr == buffer.end());
    }
  }
  // cerr << "Volume intersection totals min/max: " << voltot.min() << " " << voltot.max() << endl;

  // Now we know the total volume intersected for each of our donor nodes, and which domains need that donor info.
  // Send the donor information to each domain.
  for (auto iproc = 0; iproc != nprocs; ++iproc) {
    const unsigned n = donorNodeIDs[iproc].size();
    if (n > 0) {
      sendBuffers.push_back(vector<char>());
      for (auto i: donorNodeIDs[iproc]) {
        packElement(voltot(i), sendBuffers.back());
        for (unsigned kk = 0; kk != nScalarFields; ++kk) packElement((*scalarDonorFields[kk])(i), sendBuffers.back());
        for (unsigned kk = 0; kk != nVectorFields; ++kk) packElement((*vectorDonorFields[kk])(i), sendBuffers.back());
        for (unsigned kk = 0; kk != nTensorFields; ++kk) packElement((*tensorDonorFields[kk])(i), sendBuffers.back());
        for (unsigned kk = 0; kk != nSymTensorFields; ++kk) packElement((*symTensorDonorFields[kk])(i), sendBuffers.back());
      }
      sendBufsizes.push_back(sendBuffers.back().size());
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&sendBufsizes.back(), 1, MPI_UNSIGNED, iproc, 30, Communicator::communicator(), &sendRequests.back());
      sendRequests.push_back(MPI_Request());
      MPI_Isend(&sendBuffers.back()[0], sendBufsizes.back(), MPI_CHAR, iproc, 31, Communicator::communicator(), &sendRequests.back());
    }
  }

  // Finally get the total donor information from each domain and splat onto our acceptors.
  vector<Scalar> scalarDonorValues(nScalarFields);
  vector<Vector> vectorDonorValues(nVectorFields);
  vector<Tensor> tensorDonorValues(nTensorFields);
  vector<SymTensor> symTensorDonorValues(nSymTensorFields);

  for (auto iproc = 0; iproc != nprocs; ++iproc) {
    const unsigned ndonors = intersectDonorIndices[iproc].size();
    if (ndonors > 0) {
      unsigned bufsize;
      MPI_Status recvStatus;
      MPI_Recv(&bufsize, 1, MPI_UNSIGNED, iproc, 30, Communicator::communicator(), &recvStatus);
      CHECK(bufsize > 0);
      vector<char> buffer(bufsize);
      recvStatus = MPI_Status();
      MPI_Recv(&buffer[0], bufsize, MPI_CHAR, iproc, 31, Communicator::communicator(), &recvStatus);
      vector<char>::const_iterator bufItr = buffer.begin();
      for (unsigned i = 0; i != ndonors; ++i) {
        const unsigned nacceptors = intersectAcceptorIndices[iproc][i].size();
        double voltoti;
        unpackElement(voltoti, bufItr, buffer.end());
        for (unsigned kk = 0; kk != nScalarFields; ++kk) unpackElement(scalarDonorValues[kk], bufItr, buffer.end());
        for (unsigned kk = 0; kk != nVectorFields; ++kk) unpackElement(vectorDonorValues[kk], bufItr, buffer.end());
        for (unsigned kk = 0; kk != nTensorFields; ++kk) unpackElement(tensorDonorValues[kk], bufItr, buffer.end());
        for (unsigned kk = 0; kk != nSymTensorFields; ++kk) unpackElement(symTensorDonorValues[kk], bufItr, buffer.end());
        CHECK(voltoti > 0.0);
        for (unsigned k = 0; k != nacceptors; ++k) {
          const unsigned j = intersectAcceptorIndices[iproc][i][k];
          const double volj = intersectVols[iproc][i][k];
          const Scalar f = volj/voltoti;
          // if ((posA(j) - Vector(0.35,0.35,0.25)).magnitude() < 1e-5) {
          //   cerr << " --> " << i << " " << j << " " << voltoti << " " << volj << " " << f << endl;
          // }
          for (unsigned kk = 0; kk != nScalarFields; ++kk) (*scalarAcceptorFields[kk])(j) += f*scalarDonorValues[kk];
          for (unsigned kk = 0; kk != nVectorFields; ++kk) (*vectorAcceptorFields[kk])(j) += f*vectorDonorValues[kk];
          for (unsigned kk = 0; kk != nTensorFields; ++kk) (*tensorAcceptorFields[kk])(j) += f*tensorDonorValues[kk];
          for (unsigned kk = 0; kk != nSymTensorFields; ++kk) (*symTensorAcceptorFields[kk])(j) += f*symTensorDonorValues[kk];
        }
      }
      CHECK(bufItr == buffer.end());
    }
  }

  // Finally, wait 'til all our communication is completed before exiting.
  unsigned numSendRequests = sendRequests.size();
  vector<MPI_Status> sendStatus(numSendRequests);
  MPI_Waitall(numSendRequests, &(*sendRequests.begin()), &(*sendStatus.begin()));

#else
  //..........................................................................
  // Non-parallel version.

  // Look for intersecting node volumes.
  Field<Dimension, vector<unsigned>> intersectIndices("intersection indices", *donorNodeListPtr);
  Field<Dimension, vector<Scalar>> intersectVols("intesection volumes", *donorNodeListPtr);
  for (unsigned i = 0; i != nD; ++i) {
    vector<int> Amasters, Acoarse, Arefine;
    neighborA.setMasterList(posD(i), HD(i), Amasters, Acoarse);
    neighborA.setRefineNeighborList(posD(i), HD(i), Acoarse, Arefine);
    for (auto jitr = Arefine.begin(); jitr < Arefine.end(); ++jitr) {
      const int j = *jitr;
      if (localDonorCells(i).intersect(localAcceptorCells(j))) {
        const vector<Facet>& facets = localAcceptorCells(j).facets();
        vector<Plane> planes;
        planes.reserve(facets.size());
        for (const Facet& facet: facets) planes.push_back(Plane(facet.position(), -facet.normal()));
        const Scalar Vi = clippedVolume(localDonorCells(i), planes);
        // if (i == 451) cerr << " --> " << i << " " << j << " " << Vi << " " << localDonorCells(i).volume() << " " << localAcceptorCells(j).volume() << endl;
        if (Vi > 0.0) {
          // if (i == 0 and j == 0) cerr << "   " << i << " -> " << j << " : " << Vi << " " << Vi/localDonorCells(i).volume() << endl;
          intersectIndices(i).push_back(j);
          intersectVols(i).push_back(Vi);
        }
      }
    }
  }

  // cerr << "Node 0: " << endl
  //      << "        ";
  // copy(intersectIndices(0).begin(), intersectIndices(0).end(), std::ostream_iterator<unsigned>(cerr, " "));
  // cerr << endl
  //      << "        ";
  // copy(intersectVols(0).begin(), intersectVols(0).end(), std::ostream_iterator<double>(cerr, " "));
  // cerr << endl
  //      << "        " << std::accumulate(intersectVols(0).begin(), intersectVols(0).end(), 0.0) << endl;

  // Now we can go through and splat the conserved values from the donor to acceptor volumes.
  // double voltotmin = DBL_MAX, voltotmax = DBL_MIN;
  for (unsigned i = 0; i != nD; ++i) {
    const unsigned n = intersectIndices(i).size();
    CHECK(intersectVols(i).size() == n);
    const Scalar voltotInv = safeInv(accumulate(intersectVols(i).begin(), intersectVols(i).end(), 0.0));
    // voltotmin = min(voltotmin, accumulate(intersectVols(i).begin(), intersectVols(i).end(), 0.0));
    // voltotmax = max(voltotmax, accumulate(intersectVols(i).begin(), intersectVols(i).end(), 0.0));
    // if (accumulate(intersectVols(i).begin(), intersectVols(i).end(), 0.0) < 1.0e-4) {
    //   cerr << "Strange node: " << i << " " << posD(i) << " " << accumulate(intersectVols(i).begin(), intersectVols(i).end(), 0.0) << endl;
    // }
    for (unsigned k = 0; k != n; ++k) {
      const unsigned j = intersectIndices(i)[k];
      const Scalar f = intersectVols(i)[k]*voltotInv;
      for (unsigned kk = 0; kk != nScalarFields; ++kk) (*scalarAcceptorFields[kk])(j) += f*(*scalarDonorFields[kk])(i);
      for (unsigned kk = 0; kk != nVectorFields; ++kk) (*vectorAcceptorFields[kk])(j) += f*(*vectorDonorFields[kk])(i);
      for (unsigned kk = 0; kk != nTensorFields; ++kk) (*tensorAcceptorFields[kk])(j) += f*(*tensorDonorFields[kk])(i);
      for (unsigned kk = 0; kk != nSymTensorFields; ++kk) (*symTensorAcceptorFields[kk])(j) += f*(*symTensorDonorFields[kk])(i);
    }
  }
  // cerr << "Volume intersection totals min/max: " << voltotmin << " " << voltotmax << endl;

#endif

}

}
