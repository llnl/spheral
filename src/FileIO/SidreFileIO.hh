//---------------------------------Spheral++----------------------------------//
// SidreFileIO -- Provide the interface to sidre file objects.
//
// Created by Mikhail Zakharchanka, 11/4/2021
//----------------------------------------------------------------------------//
#ifndef __Spheral_SidreFileIO__
#define __Spheral_SidreFileIO__

#include "FileIO.hh"

#include <vector>

// Forward declarations
namespace axom {
  namespace sidre {
    class DataStore;
    class Group;
  }
}

namespace Spheral
{

class SidreFileIO: public FileIO
{
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  SidreFileIO();
  SidreFileIO(const std::string fileName, AccessType access);
  SidreFileIO(const std::string fileName, AccessType access, int numFiles);

  // Destructor.
  virtual ~SidreFileIO();

  // All File objects must provide methods to open and close the files.
  virtual void open(const std::string fileName, AccessType access) override;
  virtual void close() override;

  // Used to pass a Sidre group if you don't want to have fileIO as part of the constructor
  void setGroup(axom::sidre::Group* group);

  //******************************************************************************
  // Methods all FileIO descendent classes must provide.
  //******************************************************************************
  // Check if the specified path is in the file.
  virtual bool pathExists(const std::string path) const override;

  // All FileIO objects had better be able to read and write the primitive
  // DataTypes.
  virtual void write(const unsigned& value, const std::string path) override;
  virtual void write(const size_t& value, const std::string path) override;
  virtual void write(const int& value, const std::string path) override;
  virtual void write(const bool& value, const std::string path) override;
  virtual void write(const double& value, const std::string path) override;
  virtual void write(const std::string& value, const std::string path) override;
  virtual void write(const std::vector<int>& value, const std::string path) override;
  virtual void write(const std::vector<double>& value, const std::string path) override;
  virtual void write(const std::vector<std::string>& value, const std::string path) override;

  virtual void write(const Dim<1>::Vector& value, const std::string path) override;
  virtual void write(const Dim<1>::Tensor& value, const std::string path) override;
  virtual void write(const Dim<1>::SymTensor& value, const std::string path) override;
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string path) override;

  virtual void write(const Dim<2>::Vector& value, const std::string path) override;
  virtual void write(const Dim<2>::Tensor& value, const std::string path) override;
  virtual void write(const Dim<2>::SymTensor& value, const std::string path) override;
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string path) override;

  virtual void write(const Dim<3>::Vector& value, const std::string path) override;
  virtual void write(const Dim<3>::Tensor& value, const std::string path) override;
  virtual void write(const Dim<3>::SymTensor& value, const std::string path) override;
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string path) override;

  virtual void read(unsigned& value, const std::string path) const override;
  virtual void read(size_t& value, const std::string path) const override;
  virtual void read(int& value, const std::string path) const override;
  virtual void read(bool& value, const std::string path) const override;
  virtual void read(double& value, const std::string path) const override;
  virtual void read(std::string& value, const std::string path) const override;
  virtual void read(std::vector<int>& value, const std::string path) const override;
  virtual void read(std::vector<double>& value, const std::string path) const override;
  virtual void read(std::vector<std::string>& value, const std::string path) const override;

  virtual void read(Dim<1>::Vector& value, const std::string path) const override;
  virtual void read(Dim<1>::Tensor& value, const std::string path) const override;
  virtual void read(Dim<1>::SymTensor& value, const std::string path) const override;
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string path) const override;

  virtual void read(Dim<2>::Vector& value, const std::string path) const override;
  virtual void read(Dim<2>::Tensor& value, const std::string path) const override;
  virtual void read(Dim<2>::SymTensor& value, const std::string path) const override;
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string path) const override;

  virtual void read(Dim<3>::Vector& value, const std::string path) const override;
  virtual void read(Dim<3>::Tensor& value, const std::string path) const override;
  virtual void read(Dim<3>::SymTensor& value, const std::string path) const override;
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string path) const override;
  //******************************************************************************

  // Forward hidden base class methods
  using FileIO::write;
  using FileIO::read;

private:
  //--------------------------- Private Interface ---------------------------//
  // A pointer to the root of the sidre datastore associated with this object.
  std::unique_ptr<axom::sidre::DataStore> mDataStorePtr;

  // A pointer to the group that will have data written to it, in standalone Spheral's
  // case this will be the root group. Other codes could pass a group to Spheral to use
  // their datastore.
  axom::sidre::Group* baseGroup = nullptr;

  // write() function in sidre needs to have access to file name
  std::string mFileName;

  // The number of restart files you want SPIO (sidre parallel IO) to write
  int numRestartFiles = Process::getTotalNumberOfProcesses(); //MPI_Comm_size()

  // Don't allow assignment.
  SidreFileIO& operator=(const SidreFileIO& rhs);
};

}

#endif
