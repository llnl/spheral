include_directories(.)
set(Solvers_inst )

set(Solvers_sources
  KINSOL.cc
)

instantiate(Solvers_inst Solvers_sources)

set(Solvers_headers
  SolverFunction.hh
  KINSOL.hh
)

spheral_add_obj_library(Solvers SPHERAL_OBJ_LIBS)
