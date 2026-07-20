include_directories(.)
set(Field_inst
    AllNodeIterator
    CoarseNodeIterator
    GhostNodeIterator
    InternalNodeIterator
    MasterNodeIterator
    NodeIteratorBase
    RefineNodeIterator
   )


set(Field_sources )

instantiate(Field_inst Field_sources)

set(Field_headers
    AggregateField.hh
    AllNodeIterator.hh
    AllNodeIteratorInline.hh
    CoarseNodeIterator.hh
    CoarseNodeIteratorInline.hh
    FieldBase.hh
    FieldBaseInline.hh
    FieldView.hh
    FieldViewInline.hh
    Field.hh
    FieldInline.hh
    FieldListBase.hh
    FieldListBaseInline.hh
    FieldList.hh
    FieldListInline.hh
    FieldListView.hh
    FieldListViewInline.hh
    FieldListSet.hh
    FieldListThreadWrapper.hh
    GhostNodeIterator.hh
    GhostNodeIteratorInline.hh
    InternalNodeIterator.hh
    InternalNodeIteratorInline.hh
    MasterNodeIterator.hh
    MasterNodeIteratorInline.hh
    NodeIteratorBase.hh
    NodeIteratorBaseInline.hh
    NodeIterators.hh
    RefineNodeIterator.hh
    RefineNodeIteratorInline.hh
    uvm_allocator.hh
    )

spheral_add_obj_library(Field SPHERAL_OBJ_LIBS)

