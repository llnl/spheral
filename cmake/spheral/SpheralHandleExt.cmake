
#----------------------------------------------------------------------------------------
#                                   Spheral_Handle_Ext
#----------------------------------------------------------------------------------------

# ----------------------
# INPUT VARIABLES
# ----------------------
# <lib_name>      : REQUIRED : name of target TPL
# APPLE           : REQUIRED : flag for Mac OSX

# ----------------------
# OUTPUT VARIABLES
# ----------------------
# <lib_name>_libs : list of library names with modified extension
#----------------------------------------------------------------------------------------

function(Spheral_Handle_Ext lib_name APPLE)

  if(APPLE)
    set(SHARED_EXT "dylib")
  else()
    set(SHARED_EXT "so")
  endif()

  if(ENABLE_STATIC_TPL)
    string(REPLACE ".${SHARED_EXT}" ".a;" ${lib_name}_libs ${${lib_name}_libs})
  endif()

endfunction()
