set(Field_inst
   )

set(Field_sources
   )

instantiate(Field_inst Field_sources)

set(Field_headers
    AggregateField.hh
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
    uvm_allocator.hh
    )

spheral_add_package(Field)
