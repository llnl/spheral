
#----------------------------------------------------------------------------------------
#                                   Spheral_Handle_Ext
#----------------------------------------------------------------------------------------

# ----------------------
# INPUT VARIABLES
# ----------------------
# <lib_name>      : REQUIRED : name of target TPL

# ----------------------
# OUTPUT VARIABLES
# ----------------------
# <lib_name>_libs : list of library names with modified extension
#----------------------------------------------------------------------------------------

function(Spheral_Handle_Ext lib_name)
  if(APPLE)
    set(_lib_ext ".dylib")
  elseif(ENABLE_STATIC_TPL)
    set(_lib_ext ".a")
  else()
    list(GET ${lib_name}_libs 0 _test_lib)
    set(_full_test_lib "${${lib_name}_DIR}/**/${_test_lib}")
    file(GLOB FOUND_DYLIB "${_full_test_lib}.dylib")
    file(GLOB FOUND_STATIC "${_full_test_lib}.a")
    if(FOUND_DYLIB)
      set(_lib_ext ".dylib")
    elseif(FOUND_STATIC)
      set(_lib_ext ".a")
    else()
      set(_lib_ext ".so")
    endif()
  endif()

  list(TRANSFORM ${lib_name}_libs APPEND ${_lib_ext})
  set(${lib_name}_libs "${${lib_name}_libs}" PARENT_SCOPE)
  message("${${lib_name}_libs}")
endfunction()
