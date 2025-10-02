include_directories(.)

set(Field_sources )

if(ENABLE_INSTANTIATIONS)
  list(APPEND Field_sources
       AllNodeIteratorInst.cc
       CoarseNodeIteratorInst.cc
       GhostNodeIteratorInst.cc
       InternalNodeIteratorInst.cc
       MasterNodeIteratorInst.cc
       NodeIteratorBaseInst.cc
       RefineNodeIteratorInst.cc
       )
endif()

set(Field_headers
    AggregateField.hh
    AllNodeIterator.hh
    AllNodeIteratorInline.hh
    CoarseNodeIterator.hh
    CoarseNodeIteratorInline.hh
    Field.hh
    FieldBase.hh
    FieldBaseInline.hh
    FieldInline.hh
    FieldList.hh
    FieldListBase.hh
    FieldListBaseInline.hh
    FieldListInline.hh
    FieldListSet.hh
    FieldListSetInline.hh
    FieldListThreadWrapper.hh
    FieldListView.hh
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

# Even though these files have a .cc extension, they are template definitions
# and should be treated as headers.
list(APPEND Field_headers
     AllNodeIterator.cc
     CoarseNodeIterator.cc
     GhostNodeIterator.cc
     InternalNodeIterator.cc
     MasterNodeIterator.cc
     NodeIteratorBase.cc
     RefineNodeIterator.cc
     )

spheral_add_obj_library(Field SPHERAL_OBJ_LIBS)

