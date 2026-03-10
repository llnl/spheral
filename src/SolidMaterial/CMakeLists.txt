include_directories(.)
set(SolidMaterial_inst
    LinearPolynomialEquationOfState
    GruneisenEquationOfState
    OsborneEquationOfState
    MurnaghanEquationOfState
    TillotsonEquationOfState
    StrengthModel
    ConstantStrength
    SteinbergGuinanStrength
    JohnsonCookStrength
    CollinsStrength
    iSALEROCKStrength
    PhysicsEvolvingMaterialLibrary
    )

set(SolidMaterial_sources )

set(SolidMaterial_headers
    ANEOS.hh
    CollinsStrength.hh
    iSALEROCKStrength.hh
    ConstantStrength.hh
    Geodyn.hh
    GruneisenEquationOfState.hh
    GruneisenEquationOfStateInline.hh
    JohnsonCookStrength.hh
    LinearPolynomialEquationOfState.hh
    LinearPolynomialEquationOfStateInline.hh
    MurnaghanEquationOfState.hh
    MurnaghanEquationOfStateInline.hh
    NullStrength.hh
    OsborneEquationOfState.hh
    OsborneEquationOfStateInline.hh
    PhysicsEvolvingMaterialLibrary.hh
    PolynomialFit.hh
    PolynomialFitInline.hh
    SolidEquationOfState.hh
    SolidEquationOfStateInline.hh
    SteinbergGuinanLundStrength.hh
    SteinbergGuinanStrength.hh
    StrengthModel.hh
    TillotsonEquationOfState.hh
    TillotsonEquationOfStateInline.hh
    )

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/ShadowANEOS.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/ShadowANEOS.py)

spheral_install_python_files(
    MaterialPropertiesLib.py
    SolidMaterialUnits.py
    SolidMaterialEquationsOfState.py
    ShadowTillotsonEquationOfState.py
    ShadowGruneisenEquationOfState.py
    ShadowConstantStrength.py
    ${CMAKE_CURRENT_BINARY_DIR}/ShadowANEOS.py
    )

# ANEOS
if (SPHERAL_ENABLE_ANEOS)
  message("-- ANEOS enabled")
  list(APPEND SolidMaterial_inst      ANEOS)
  list(APPEND SolidMaterial_headers   ANEOS.hh)
  list(APPEND SolidMaterial_sources   ANEOS_initialize.f)
else()
  message("-- ANEOS disabled")
endif()

instantiate(SolidMaterial_inst SolidMaterial_sources)

# Ignore -W-maybe-uninitialized warnings thrown from bisectSearch for possibly uninitialized iterators from boost multi-array.
set_source_files_properties(ANEOSInst1d.cc ANEOSInst2d.cc ANEOSInst3d.cc PROPERTIES COMPILE_FLAGS "$<$<COMPILE_LANGUAGE:CXX>:-Wno-uninitialized>")

spheral_add_obj_library(SolidMaterial SPHERAL_OBJ_LIBS)

