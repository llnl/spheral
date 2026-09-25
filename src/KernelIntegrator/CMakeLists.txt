include_directories(.)
set(KernelIntegrator_inst
  BilinearIndex
  FlatConnectivity
  KernelIntegral
  KernelIntegrator
  ManufacturedSolution
  RKIntegrationKernel
  SPHIntegrationKernel
  )

set(KernelIntegrator_sources
  GaussLegendreValues.cc
  SymmetricTetrahedralValues.cc
  SymmetricTriangularValues.cc
  )
instantiate(KernelIntegrator_inst KernelIntegrator_sources)

set(KernelIntegrator_headers
  BilinearIndex.hh
  BilinearIndexInline.hh
  DataTypeFunctions.hh
  FlatConnectivity.hh
  FlatConnectivityInline.hh
  GaussLegendreValues.hh
  IntegrationCoefficient.hh
  IntegrationCoefficientInline.hh
  IntegrationKernel.hh
  KernelIntegral.hh
  KernelIntegralInline.hh
  KernelIntegrationData.hh
  KernelIntegrator.hh
  KernelIntegratorInline.hh
  ManufacturedSolution.hh
  RKIntegrationKernel.hh
  RKIntegrationKernelInline.hh
  SPHIntegrationKernel.hh
  SymmetricTetrahedralValues.hh
  SymmetricTriangularValues.hh
  )

spheral_add_obj_library(KernelIntegrator)
