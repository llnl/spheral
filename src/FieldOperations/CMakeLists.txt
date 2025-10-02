include_directories(.)

set(FieldOperations_sources )

if(ENABLE_INSTANTIATIONS)
  list(APPEND FieldOperations_sources
       smoothFieldsInst.cc
       gradientInst.cc
       divergenceInst.cc
       smoothFieldsMashInst.cc
       smoothFieldsMash2Inst.cc
       sampleFieldsMashInst.cc
       splatFieldsMashInst.cc
       gradientMashInst.cc
       divergenceMashInst.cc
       gradientPairWiseInst.cc
       divergencePairWiseInst.cc
       sampleMultipleFieldsMashInst.cc
       splatMultipleFieldsMashInst.cc
       gradDivVectorFieldListInst.cc
       gradDivVectorFieldListSimpleInst.cc
       gradDivVectorFieldListGoldenInst.cc
       gradDivVectorFieldListGolden2Inst.cc
       gradDivVectorFieldListMashInst.cc
       gradDivVectorFieldListPairWiseInst.cc
       sampleMultipleFields2LatticeInst.cc
       sampleMultipleFields2LatticeMashInst.cc
       binFieldList2LatticeInst.cc
       )
endif()

set(FieldOperations_headers
    FieldListFunctions.hh
    FieldListFunctionsMash.hh
    FieldListSecondDerivatives.hh
    PairWiseFieldListFunctions.hh
    binFieldList2Lattice.hh
    incrementPairWiseDelta.hh
    monotonicallyLimitedGradient.hh
    sampleMultipleFields2Lattice.hh
    )

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND FieldOperations_headers
     smoothFields.cc
     gradient.cc
     divergence.cc
     smoothFieldsMash.cc
     smoothFieldsMash2.cc
     sampleFieldsMash.cc
     splatFieldsMash.cc
     gradientMash.cc
     divergenceMash.cc
     gradientPairWise.cc
     divergencePairWise.cc
     sampleMultipleFieldsMash.cc
     splatMultipleFieldsMash.cc
     gradDivVectorFieldList.cc
     gradDivVectorFieldListSimple.cc
     gradDivVectorFieldListGolden.cc
     gradDivVectorFieldListGolden2.cc
     gradDivVectorFieldListMash.cc
     gradDivVectorFieldListPairWise.cc
     sampleMultipleFields2Lattice.cc
     sampleMultipleFields2LatticeMash.cc
     binFieldList2Lattice.cc
     )

spheral_add_obj_library(FieldOperations SPHERAL_OBJ_LIBS)

