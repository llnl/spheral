include_directories(.)
set(VoronoiCells_inst
  VoronoiCells
  IncrementVoronoiCells
  UpdateVoronoiCells
  SubPointPressureHourglassControl
)

set(VoronoiCells_sources
  computeVoronoiVolume.cc
  computeVoronoiVolume1d.cc
)

instantiate(VoronoiCells_inst VoronoiCells_sources)

set(VoronoiCells_headers
  VoronoiCells.hh
  computeVoronoiVolume.hh
  SubPointPressureHourglassControl.hh
)

spheral_add_obj_library(VoronoiCells SPHERAL_OBJ_LIBS)
