include_directories(.)
set(NodeList_inst
    FluidNodeList
    NodeList
    NodeListRegistrar
    SolidNodeList
    generateVoidNodes
    nthNodalMoment
    DEMNodeList
   )

set(NodeList_sources )

instantiate(NodeList_inst NodeList_sources)

set(NodeList_headers
    FluidNodeList.hh
    FluidNodeListInline.hh
    FluidNodeTraits.hh
    NodeList.hh
    NodeListRegistrar.hh
    NodeListRegistrarInline.hh
    SolidNodeList.hh
    SolidNodeListInline.hh
    generateVoidNodes.hh
    nthNodalMoment.hh
    DEMNodeList.hh
    DEMNodeListInline.hh
    )

spheral_install_python_files(
  FluidNodeLists.py
  SolidNodeLists.py
  VoidNodeLists.py
  DEMNodeLists.py
  )

spheral_add_obj_library(NodeList SPHERAL_OBJ_LIBS)

