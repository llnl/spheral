include_directories(.)
set(FieldOperations_inst
    smoothFields
    gradient
    divergence
    smoothFieldsMash
    smoothFieldsMash2
    sampleFieldsMash
    splatFieldsMash
    gradientMash
    divergenceMash
    gradientPairWise
    divergencePairWise
    sampleMultipleFieldsMash
    splatMultipleFieldsMash
    gradDivVectorFieldList
    gradDivVectorFieldListSimple
    gradDivVectorFieldListGolden
    gradDivVectorFieldListGolden2
    gradDivVectorFieldListMash
    gradDivVectorFieldListPairWise
    sampleMultipleFields2Lattice
    sampleMultipleFields2LatticeMash
    binFieldList2Lattice
   )


set(FieldOperations_sources )

instantiate(FieldOperations_inst FieldOperations_sources)

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

spheral_add_obj_library(FieldOperations)

