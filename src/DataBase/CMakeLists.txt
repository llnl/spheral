include_directories(.)

set(DataBase_sources )

if(ENABLE_INSTANTIATIONS)
  list(APPEND DataBase_sources
       DataBaseInst.cc
       ReplaceWithRatioPolicyInst.cc
       StateInst.cc
       StateBaseInst.cc
       StateDerivativesInst.cc
       )
endif()

set(DataBase_headers
    DataBase.hh
    DataBaseInline.hh
    UpdatePolicyBase.hh
    UpdatePolicyBaseInline.hh
    FieldUpdatePolicy.hh
    FieldUpdatePolicyInline.hh
    CopyState.hh
    CopyStateInline.hh
    make_policy.hh
    IncrementState.hh
    IncrementStateInline.hh
    IncrementBoundedState.hh
    IncrementBoundedStateInline.hh
    MaxReplaceState.hh
    PureReplaceState.hh
    PureReplaceStateInline.hh
    ReplaceState.hh
    ReplaceStateInline.hh
    PureReplaceBoundedState.hh
    PureReplaceBoundedStateInline.hh
    ReplaceBoundedState.hh
    ReplaceBoundedStateInline.hh
    ReplaceWithRatioPolicy.hh
    State.hh
    StateBase.hh
    StateBaseInline.hh
    StateDerivatives.hh
    StateInline.hh
    applyPolicyToFieldList.hh
    updateStateFields.hh
    )

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND DataBase_headers
     DataBase.cc
     ReplaceWithRatioPolicy.cc
     State.cc
     StateBase.cc
     StateDerivatives.cc
     )

spheral_add_obj_library(DataBase SPHERAL_OBJ_LIBS)

