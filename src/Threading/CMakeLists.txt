include_directories(.)
set(Threading_inst
  )

set(Threading_sources
    GPUUtils.cc
  )

instantiate(Threading_inst Threading_sources)

set(Threading_headers
    GPUUtils.hh
    OpenMP_wrapper.hh
    ViewManager.hh
    ViewManagerInline.hh
    ManagedViewBase.hh
    ManagedView.hh
  )

spheral_add_obj_library(Threading SPHERAL_OBJ_LIBS)

