include_directories(.)
set(Material_inst
    EquationOfState
    GammaLawGas
    PolytropicEquationOfState
    IsothermalEquationOfState
    StiffenedGas
   )

set(Material_sources
    PhysicalConstants.cc
    UnitConversion.cc
    )

set(Material_headers
    EquationOfState.hh
    EquationOfStateInline.hh
    GammaLawGas.hh
    StiffenedGas.hh
    HelmholtzEquationOfState.hh
    IsothermalEquationOfState.hh
    IsothermalEquationOfStateInline.hh
    PhysicalConstants.hh
    PhysicalConstantsInline.hh
    PolytropicEquationOfState.hh
    PolytropicEquationOfStateInline.hh
    UnitConversion.hh
    )

spheral_install_python_files(
  MaterialUnits.py
  MaterialEquationsOfState.py
  computeHugoniot.py
  )

# Helmholtz
if (SPHERAL_ENABLE_HELMHOLTZ)
  message("-- Helmholtz equation of state enabled")
  list(APPEND Material_inst      HelmholtzEquationOfState)
  list(APPEND Material_headers   HelmholtzEquationOfState.hh)
  list(APPEND Material_sources   public_helm.f90 invert_helm.f90)
else()
  message("-- Helmholtz equation of state disabled")
endif()

instantiate(Material_inst Material_sources)

spheral_add_obj_library(Material SPHERAL_OBJ_LIBS)

