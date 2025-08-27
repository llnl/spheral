include_directories(.)
set(Solvers_inst )

set(Solvers_sources
  EigenLinearSolver.cc
  HypreFunctions.cc
  HypreLinearSolver.cc
  KINSOL.cc
  LinearSolver.cc
)

instantiate(Solvers_inst Solvers_sources)

set(Solvers_headers
  containsConstantNodes.hh
  containsConstantNodesInline.hh
  EigenLinearSolver.hh
  EigenOptions.hh
  HypreFunctions.hh
  HypreLinearSolver.hh
  HypreLinearSolverInline.hh
  HypreOptions.hh
  IncrementalStatistic.hh
  IncrementalStatisticInline.hh
  KINSOL.hh
  LinearSolver.hh
  LinearSolverInline.hh
  SolverFunction.hh
)

spheral_add_obj_library(Solvers SPHERAL_OBJ_LIBS)
