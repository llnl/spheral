//---------------------------------Spheral++----------------------------------//
// RegisterMPIDataTypes
// A singleton helper class that registers special Spheral defined data types
// with MPI.
//
// Created by J. Michael Owen, Tue Jul  7 14:02:12 PDT 2009
//----------------------------------------------------------------------------//

#include "RegisterMPIDataTypes.hh"
#include "Utilities/DataTypeTraits.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor (private).
//------------------------------------------------------------------------------
RegisterMPIDataTypes::
RegisterMPIDataTypes() {

#ifdef SPHERAL_ENABLE_MPI
  // Vectors.
  MPI_Type_contiguous(DataTypeTraits<Dim<1>::Vector>::numElements(Dim<1>::Vector::zero), MPI_DOUBLE, &MPI_Vector1d);
  MPI_Type_contiguous(DataTypeTraits<Dim<2>::Vector>::numElements(Dim<2>::Vector::zero), MPI_DOUBLE, &MPI_Vector2d);
  MPI_Type_contiguous(DataTypeTraits<Dim<3>::Vector>::numElements(Dim<3>::Vector::zero), MPI_DOUBLE, &MPI_Vector3d);
  MPI_Type_commit(&MPI_Vector1d);
  MPI_Type_commit(&MPI_Vector2d);
  MPI_Type_commit(&MPI_Vector3d);

//   {
//     int block_lengths[2];
//     MPI_Aint displacements[2];
//     MPI_Aint addresses[3];
//     MPI_Datatype type_list[2];

//     // First specify the types.
//     type_list[0] = MPI_DOUBLE;
//     type_list[1] = MPI_DOUBLE;

//     // Specify the number of elements of each type.
//     block_lengths[0] = 1;
//     block_lengths[1] = 1;

//     // Calculate the displacements of the members
//     // relative to indata.
//     Dim<2>::Vector tmp;
//     MPI_Address(&tmp, &addresses[0]);
//     MPI_Address(&tmp(0), &addresses[1]);
//     MPI_Address(&tmp(1), &addresses[2]);
//     displacements[0] = addresses[1] - addresses[0];
//     displacements[1] = addresses[2] - addresses[0];

//     // Create the derived type.
//     MPI_Type_struct(2, block_lengths, displacements, type_list, &MPI_Vector2d);
//     MPI_Type_commit(&MPI_Vector2d);
//   }

  // Tensors.
  MPI_Type_contiguous(DataTypeTraits<Dim<1>::Tensor>::numElements(Dim<1>::Tensor::zero), MPI_DOUBLE, &MPI_Tensor1d);
  MPI_Type_contiguous(DataTypeTraits<Dim<2>::Tensor>::numElements(Dim<2>::Tensor::zero), MPI_DOUBLE, &MPI_Tensor2d);
  MPI_Type_contiguous(DataTypeTraits<Dim<3>::Tensor>::numElements(Dim<3>::Tensor::zero), MPI_DOUBLE, &MPI_Tensor3d);
  MPI_Type_commit(&MPI_Tensor1d);
  MPI_Type_commit(&MPI_Tensor2d);
  MPI_Type_commit(&MPI_Tensor3d);

  // SymTensors.
  MPI_Type_contiguous(DataTypeTraits<Dim<1>::SymTensor>::numElements(Dim<1>::SymTensor::zero), MPI_DOUBLE, &MPI_SymTensor1d);
  MPI_Type_contiguous(DataTypeTraits<Dim<2>::SymTensor>::numElements(Dim<2>::SymTensor::zero), MPI_DOUBLE, &MPI_SymTensor2d);
  MPI_Type_contiguous(DataTypeTraits<Dim<3>::SymTensor>::numElements(Dim<3>::SymTensor::zero), MPI_DOUBLE, &MPI_SymTensor3d);
  MPI_Type_commit(&MPI_SymTensor1d);
  MPI_Type_commit(&MPI_SymTensor2d);
  MPI_Type_commit(&MPI_SymTensor3d);

  // ThirdRankTensors.
  MPI_Type_contiguous(DataTypeTraits<Dim<1>::ThirdRankTensor>::numElements(Dim<1>::ThirdRankTensor::zero), MPI_DOUBLE, &MPI_ThirdRankTensor1d);
  MPI_Type_contiguous(DataTypeTraits<Dim<2>::ThirdRankTensor>::numElements(Dim<2>::ThirdRankTensor::zero), MPI_DOUBLE, &MPI_ThirdRankTensor2d);
  MPI_Type_contiguous(DataTypeTraits<Dim<3>::ThirdRankTensor>::numElements(Dim<3>::ThirdRankTensor::zero), MPI_DOUBLE, &MPI_ThirdRankTensor3d);
  MPI_Type_commit(&MPI_ThirdRankTensor1d);
  MPI_Type_commit(&MPI_ThirdRankTensor2d);
  MPI_Type_commit(&MPI_ThirdRankTensor3d);

  // FourthRankTensors.
  MPI_Type_contiguous(DataTypeTraits<Dim<1>::FourthRankTensor>::numElements(Dim<1>::FourthRankTensor::zero), MPI_DOUBLE, &MPI_FourthRankTensor1d);
  MPI_Type_contiguous(DataTypeTraits<Dim<2>::FourthRankTensor>::numElements(Dim<2>::FourthRankTensor::zero), MPI_DOUBLE, &MPI_FourthRankTensor2d);
  MPI_Type_contiguous(DataTypeTraits<Dim<3>::FourthRankTensor>::numElements(Dim<3>::FourthRankTensor::zero), MPI_DOUBLE, &MPI_FourthRankTensor3d);
  MPI_Type_commit(&MPI_FourthRankTensor1d);
  MPI_Type_commit(&MPI_FourthRankTensor2d);
  MPI_Type_commit(&MPI_FourthRankTensor3d);

  // FifthRankTensors.
  MPI_Type_contiguous(DataTypeTraits<Dim<1>::FifthRankTensor>::numElements(Dim<1>::FifthRankTensor::zero), MPI_DOUBLE, &MPI_FifthRankTensor1d);
  MPI_Type_contiguous(DataTypeTraits<Dim<2>::FifthRankTensor>::numElements(Dim<2>::FifthRankTensor::zero), MPI_DOUBLE, &MPI_FifthRankTensor2d);
  MPI_Type_contiguous(DataTypeTraits<Dim<3>::FifthRankTensor>::numElements(Dim<3>::FifthRankTensor::zero), MPI_DOUBLE, &MPI_FifthRankTensor3d);
  MPI_Type_commit(&MPI_FifthRankTensor1d);
  MPI_Type_commit(&MPI_FifthRankTensor2d);
  MPI_Type_commit(&MPI_FifthRankTensor3d);
#endif

}

//------------------------------------------------------------------------------
// Destructor (private)
//------------------------------------------------------------------------------
RegisterMPIDataTypes::
~RegisterMPIDataTypes() {
}

}
