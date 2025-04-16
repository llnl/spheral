include_directories(.)
set(Integrator_inst
  Integrator
  ImplicitIntegrator
  PredictorCorrector 
  ForwardEuler
  SynchronousRK2
  SynchronousRK4
  CheapSynchronousRK2
  Verlet
  CrankNicolson
)

set(Integrator_sources )

if (ENABLE_SUNDIALS)
  list(APPEND Integrator_inst
    BackwardEuler
    ImplicitIntegrationVectorOperator
  )
endif()

instantiate(Integrator_inst Integrator_sources)

set(Integrator_headers
  Integrator.hh
  ImplicitIntegrator.hh
  CheapSynchronousRK2.hh
  PredictorCorrector.hh
  ForwardEuler.hh
  SynchronousRK2.hh
  SynchronousRK4.hh
  Verlet.hh
  BackwardEuler.hh
  CrankNicolson.hh
  ImplicitIntegrationVectorOperator.hh
)

spheral_add_obj_library(Integrator SPHERAL_OBJ_LIBS)

