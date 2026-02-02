//---------------------------------Spheral++----------------------------------//
// PyFileIO -- A python friendly version of the FileIO interface, for use
// creating python FileIO objects.
//
// This class overrides the FileIO::read methods, since BPL has problems 
// handing back references to int, bool, and double.
//
// Created by JMO, Tue Dec 27 22:12:00 PST 2005
//----------------------------------------------------------------------------//
#ifndef __Spheral_PyFileIO_hh__
#define __Spheral_PyFileIO_hh__

#include "FileIO/FileIO.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/DBC.hh"

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include <vector>

namespace Spheral {

namespace py = pybind11;

class PyFileIO: public FileIO {
public:
  //--------------------------- Public Interface ---------------------------//
  // Constructors.
  PyFileIO();
  PyFileIO(const std::string filename, AccessType access);

  // Destructor.
  virtual ~PyFileIO();

  // Check if the specified path is in the file.
  virtual bool pathExists(const std::string) const override { VERIFY2(false, "pathExists not overridden"); return false; }

  //**************************************************************************
  // Descendent python objects can optionally override these disambiguated methods
  // Write
  virtual void write_unsigned_int(const unsigned value, const std::string path)                  override { this->writeObject(py::cast(value), path); }
  virtual void write_size_t(const size_t value, const std::string path)                          override { this->writeObject(py::cast(value), path); }
  virtual void write_int(const int value, const std::string path)                                override { this->writeObject(py::cast(value), path); }
  virtual void write_bool(const bool value, const std::string path)                              override { this->writeObject(py::cast(value), path); }
  virtual void write_double(const double value, const std::string path)                          override { this->writeObject(py::cast(value), path); }
  virtual void write_string(const std::string value, const std::string path)                     override { this->writeObject(py::cast(value), path); }
  virtual void write_vector_char(const std::vector<char>& value, const std::string path)         override;
  virtual void write_vector_int(const std::vector<int>& value, const std::string path)                    { this->writeObject(py::cast(value), path); }
  virtual void write_vector_double(const std::vector<double>& value, const std::string path)              { this->writeObject(py::cast(value), path); }
  virtual void write_vector_string(const std::vector<std::string>& value, const std::string path)         { this->writeObject(py::cast(value), path); }

  virtual void write_Vector1d(const Dim<1>::Vector& value, const std::string path)                        { this->writeObject(py::cast(value), path); }
  virtual void write_Tensor1d(const Dim<1>::Tensor& value, const std::string path)                        { this->writeObject(py::cast(value), path); }
  virtual void write_SymTensor1d(const Dim<1>::SymTensor& value, const std::string path)                  { this->writeObject(py::cast(value), path); }
  virtual void write_ThirdRankTensor1d(const Dim<1>::ThirdRankTensor& value, const std::string path)      { this->writeObject(py::cast(value), path); }
                                                                                                                                       
  virtual void write_Vector2d(const Dim<2>::Vector& value, const std::string path)                        { this->writeObject(py::cast(value), path); }
  virtual void write_Tensor2d(const Dim<2>::Tensor& value, const std::string path)                        { this->writeObject(py::cast(value), path); }
  virtual void write_SymTensor2d(const Dim<2>::SymTensor& value, const std::string path)                  { this->writeObject(py::cast(value), path); }
  virtual void write_ThirdRankTensor2d(const Dim<2>::ThirdRankTensor& value, const std::string path)      { this->writeObject(py::cast(value), path); }
                                                                                                                                       
  virtual void write_Vector3d(const Dim<3>::Vector& value, const std::string path)                        { this->writeObject(py::cast(value), path); }
  virtual void write_Tensor3d(const Dim<3>::Tensor& value, const std::string path)                        { this->writeObject(py::cast(value), path); }
  virtual void write_SymTensor3d(const Dim<3>::SymTensor& value, const std::string path)                  { this->writeObject(py::cast(value), path); }
  virtual void write_ThirdRankTensor3d(const Dim<3>::ThirdRankTensor& value, const std::string path)      { this->writeObject(py::cast(value), path); }

  // Read
  virtual unsigned read_unsigned_int(const std::string path)                               const override { return this->readObject(path).cast<unsigned>(); }                       
  virtual size_t read_size_t(const std::string path)                                       const override { return this->readObject(path).cast<size_t>(); }     
  virtual int read_int(const std::string path)                                             const override { return this->readObject(path).cast<int>(); }        
  virtual bool read_bool(const std::string path)                                           const override { return this->readObject(path).cast<bool>(); }       
  virtual double read_double(const std::string path)                                       const override { return this->readObject(path).cast<double>(); }     
  virtual std::string read_string(const std::string path)                                  const override { return this->readObject(path).cast<std::string>(); }
  virtual std::vector<char> read_vector_char(const std::string path)                       const override;
  virtual void read_vector_int(std::vector<int>& value, const std::string path)            const          { value = this->readObject(path).cast<std::vector<int>>(); }
  virtual void read_vector_double(std::vector<double>& value, const std::string path)      const          { value = this->readObject(path).cast<std::vector<double>>(); }
  virtual void read_vector_string(std::vector<std::string>& value, const std::string path) const          { value = this->readObject(path).cast<std::vector<std::string>>(); }

  virtual Dim<1>::Vector read_Vector1d(const std::string path)                             const          { return this->readObject(path).cast<Dim<1>::Vector>(); }         
  virtual Dim<1>::Tensor read_Tensor1d(const std::string path)                             const          { return this->readObject(path).cast<Dim<1>::Tensor>(); }         
  virtual Dim<1>::SymTensor read_SymTensor1d(const std::string path)                       const          { return this->readObject(path).cast<Dim<1>::SymTensor>(); }      
  virtual Dim<1>::ThirdRankTensor read_ThirdRankTensor1d(const std::string path)           const          { return this->readObject(path).cast<Dim<1>::ThirdRankTensor>(); }
                                                                                                                                                                            
  virtual Dim<2>::Vector read_Vector2d(const std::string path)                             const          { return this->readObject(path).cast<Dim<2>::Vector>(); }         
  virtual Dim<2>::Tensor read_Tensor2d(const std::string path)                             const          { return this->readObject(path).cast<Dim<2>::Tensor>(); }         
  virtual Dim<2>::SymTensor read_SymTensor2d(const std::string path)                       const          { return this->readObject(path).cast<Dim<2>::SymTensor>(); }      
  virtual Dim<2>::ThirdRankTensor read_ThirdRankTensor2d(const std::string path)           const          { return this->readObject(path).cast<Dim<2>::ThirdRankTensor>(); }
                                                                                                                                                                            
  virtual Dim<3>::Vector read_Vector3d(const std::string path)                             const          { return this->readObject(path).cast<Dim<3>::Vector>(); }         
  virtual Dim<3>::Tensor read_Tensor3d(const std::string path)                             const          { return this->readObject(path).cast<Dim<3>::Tensor>(); }         
  virtual Dim<3>::SymTensor read_SymTensor3d(const std::string path)                       const          { return this->readObject(path).cast<Dim<3>::SymTensor>(); }      
  virtual Dim<3>::ThirdRankTensor read_ThirdRankTensor3d(const std::string path)           const          { return this->readObject(path).cast<Dim<3>::ThirdRankTensor>(); }

  //***************************************************************************
  // Override the base FileIO read methods to use the above virtual methods.
  // Write methods.
  virtual void write(const unsigned& value, const std::string path)                              override { write_unsigned_int(value, path); } 
  virtual void write(const size_t& value, const std::string path)                                override { write_size_t(value, path); }       
  virtual void write(const int& value, const std::string path)                                   override { write_int(value, path); }          
  virtual void write(const bool& value, const std::string path)                                  override { write_bool(value, path); }         
  virtual void write(const double& value, const std::string path)                                override { write_double(value, path); }       
  virtual void write(const std::string& value, const std::string path)                           override { write_string(value, path); }       
  virtual void write(const std::vector<int>& value, const std::string path)                      override { write_vector_int(value, path); }   
  virtual void write(const std::vector<double>& value, const std::string path)                   override { write_vector_double(value, path); }
  virtual void write(const std::vector<std::string>& value, const std::string path)              override { write_vector_string(value, path); }

  virtual void write(const Dim<1>::Vector& value, const std::string path)                        override { write_Vector1d(value, path); }
  virtual void write(const Dim<1>::Tensor& value, const std::string path)                        override { write_Tensor1d(value, path); }
  virtual void write(const Dim<1>::SymTensor& value, const std::string path)                     override { write_SymTensor1d(value, path); }
  virtual void write(const Dim<1>::ThirdRankTensor& value, const std::string path)               override { write_ThirdRankTensor1d(value, path); }

  virtual void write(const Dim<2>::Vector& value, const std::string path)                        override { write_Vector2d(value, path); }
  virtual void write(const Dim<2>::Tensor& value, const std::string path)                        override { write_Tensor2d(value, path); }
  virtual void write(const Dim<2>::SymTensor& value, const std::string path)                     override { write_SymTensor2d(value, path); }
  virtual void write(const Dim<2>::ThirdRankTensor& value, const std::string path)               override { write_ThirdRankTensor2d(value, path); }

  virtual void write(const Dim<3>::Vector& value, const std::string path)                        override { write_Vector3d(value, path); }
  virtual void write(const Dim<3>::Tensor& value, const std::string path)                        override { write_Tensor3d(value, path); }
  virtual void write(const Dim<3>::SymTensor& value, const std::string path)                     override { write_SymTensor3d(value, path); }
  virtual void write(const Dim<3>::ThirdRankTensor& value, const std::string path)               override { write_ThirdRankTensor3d(value, path); }

  // Read methods.
  virtual void read(unsigned& value, const std::string path)                               const override { value = read_unsigned_int(path); }
  virtual void read(size_t& value, const std::string path)                                 const override { value = read_size_t(path); }       
  virtual void read(int& value, const std::string path)                                    const override { value = read_int(path); }          
  virtual void read(bool& value, const std::string path)                                   const override { value = read_bool(path); }         
  virtual void read(double& value, const std::string path)                                 const override { value = read_double(path); }       
  virtual void read(std::string& value, const std::string path)                            const override { value = read_string(path); }       
  virtual void read(std::vector<int>& value, const std::string path)                       const override { read_vector_int(value, path); }   
  virtual void read(std::vector<double>& value, const std::string path)                    const override { read_vector_double(value, path); }
  virtual void read(std::vector<std::string>& value, const std::string path)               const override { read_vector_string(value, path); }

  virtual void read(Dim<1>::Vector& value, const std::string path)                         const override { value = read_Vector1d(path); }
  virtual void read(Dim<1>::Tensor& value, const std::string path)                         const override { value = read_Tensor1d(path); }
  virtual void read(Dim<1>::SymTensor& value, const std::string path)                      const override { value = read_SymTensor1d(path); }
  virtual void read(Dim<1>::ThirdRankTensor& value, const std::string path)                const override { value = read_ThirdRankTensor1d(path); }

  virtual void read(Dim<2>::Vector& value, const std::string path)                         const override { value = read_Vector2d(path); }
  virtual void read(Dim<2>::Tensor& value, const std::string path)                         const override { value = read_Tensor2d(path); }
  virtual void read(Dim<2>::SymTensor& value, const std::string path)                      const override { value = read_SymTensor2d(path); }
  virtual void read(Dim<2>::ThirdRankTensor& value, const std::string path)                const override { value = read_ThirdRankTensor2d(path); }

  virtual void read(Dim<3>::Vector& value, const std::string path)                         const override { value = read_Vector3d(path); }
  virtual void read(Dim<3>::Tensor& value, const std::string path)                         const override { value = read_Tensor3d(path); }
  virtual void read(Dim<3>::SymTensor& value, const std::string path)                      const override { value = read_SymTensor3d(path); }
  virtual void read(Dim<3>::ThirdRankTensor& value, const std::string path)                const override { value = read_ThirdRankTensor3d(path); }
  //***************************************************************************

  // //------------------------------------------------------------------------------
  // // We have to forward the templated write/read methods to the base class due to
  // // function hiding.
  // // Write/read a vector<Value> if Value is a primitive we already know about.
  // template<typename Value> void write(const std::vector<Value>& x, const std::string path) { FileIO::write(x, path); }
  // template<typename Value> void  read(std::vector<Value>& x, const std::string path) const { FileIO::read(x, path); }
  // //------------------------------------------------------------------------------

  // Forward hidden base class methods
  using FileIO::write;
  using FileIO::read;

private:
  //--------------------------- Private Interface ---------------------------//
  // Don't allow assignment.
  PyFileIO& operator=(const FileIO& rhs);
};

}

#endif
