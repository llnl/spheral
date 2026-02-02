//---------------------------------Spheral++----------------------------------//
// siloMeshDump
//
// A set of methods for writing out Spheral++ data in meshed format to silo
// files.
//
// Created by JMO, Wed Dec  8 22:09:19 PST 2010
//----------------------------------------------------------------------------//
#include <algorithm>
#include <set>

#include "silo.h"

#include "Mesh.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Write a master silo file.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
writeSiloMeshMasterFile(const string baseName,
                        const Mesh<Dimension>& mesh,
                        const string label,
                        const unsigned numDomains) {

  // Create the silo file.
  string fileName = baseName + ".silo";
  DBfile* db = DBCreate(fileName.c_str(), DB_CLOBBER, DB_LOCAL, label.c_str(), DB_HDF5);

  // Write the domain file list.
  char* mmeshName = "MMESH";
  vector<char*> meshNames;
  for (unsigned idomain = 0; idomain != numDomains; ++idomain) 
    meshNames.push_back(
  VERIFY(DBPutMultiMesh(db, (char*) "MMESH", 

  // Clean up.
  db.Close();
  delete db;
}

//------------------------------------------------------------------------------
// The top method the user should call.
//------------------------------------------------------------------------------
template<typename Dimension>
void
siloMeshDump(const string baseName,
             const Mesh<Dimension>& mesh,
             const string label) {

  // Parallel info.
  const unsigned domainID = Process::getRank();
  const unsigned numDomains = Process::getTotalNumberOfProcesses();

  // If we're process 0 write the master file.
  if (domainID == 0) writeSiloMeshMasterFile(baseName, mesh, label, numDomains);

  // All processes write their own info.
  writeSiloDomainFile(baseName, domainID, numDomains);

  // That's it.
}

}
