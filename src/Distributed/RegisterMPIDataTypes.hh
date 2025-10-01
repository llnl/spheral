//---------------------------------Spheral++----------------------------------//
// RegisterMPIDataTypes
// A singleton helper class that registers special Spheral defined data types
// with MPI.
//
// Created by J. Michael Owen, Tue Jul  7 14:02:12 PDT 2009
//----------------------------------------------------------------------------//
#ifndef __Spheral__RegisterMPIDataTypes__
#define __Spheral__RegisterMPIDataTypes__

#include "Geometry/Dimension.hh"
#include "config.hh"
#ifdef ENABLE_MPI
#include <mpi.h>
#endif

namespace Spheral {

class RegisterMPIDataTypes {
public:

  // Get the instance.
  static RegisterMPIDataTypes& instance();

  // The MPI types we're here to provide.
#ifdef ENABLE_MPI
  MPI_Datatype MPI_Vector1d, MPI_Vector2d, MPI_Vector3d;
  MPI_Datatype MPI_Tensor1d, MPI_Tensor2d, MPI_Tensor3d;
  MPI_Datatype MPI_SymTensor1d, MPI_SymTensor2d, MPI_SymTensor3d;
  MPI_Datatype MPI_ThirdRankTensor1d, MPI_ThirdRankTensor2d, MPI_ThirdRankTensor3d;
  MPI_Datatype MPI_FourthRankTensor1d, MPI_FourthRankTensor2d, MPI_FourthRankTensor3d;
  MPI_Datatype MPI_FifthRankTensor1d, MPI_FifthRankTensor2d, MPI_FifthRankTensor3d;
#endif

private:
  RegisterMPIDataTypes();
  RegisterMPIDataTypes(const RegisterMPIDataTypes&);
  RegisterMPIDataTypes& operator=(const RegisterMPIDataTypes&);
  ~RegisterMPIDataTypes();
};

}

#include "RegisterMPIDataTypesInline.hh"

#endif

