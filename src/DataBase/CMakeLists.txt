include_directories(.)
set(DataBase_inst
    DataBase
    ReplaceWithRatioPolicy
    State
    StateBase
    StateDerivatives
   )


set(DataBase_sources )

instantiate(DataBase_inst DataBase_sources)

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

spheral_add_obj_library(DataBase SPHERAL_OBJ_LIBS)

