if(NOT dune-distrib-test_FOUND)

####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was dune-distrib-test-config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

#report other information
set_and_check(dune-distrib-test_PREFIX "${PACKAGE_PREFIX_DIR}")
set_and_check(dune-distrib-test_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")
set(dune-distrib-test_CXX_FLAGS " -std=c++0x ")
set(dune-distrib-test_CXX_FLAGS_DEBUG "-g -std=c++0x ")
set(dune-distrib-test_CXX_FLAGS_MINSIZEREL "-Os -DNDEBUG -std=c++0x ")
set(dune-distrib-test_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -std=c++0x ")
set(dune-distrib-test_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG -std=c++0x ")
set(dune-distrib-test_DEPENDS "dune-grid;dune-pdelab")
set(dune-distrib-test_SUGGESTS "")
set(dune-distrib-test_MODULE_PATH "${PACKAGE_PREFIX_DIR}/share/dune/cmake/modules")
set(dune-distrib-test_LIBRARIES "")
#import the target
if(dune-distrib-test_LIBRARIES)
  get_filename_component(_dir "${CMAKE_CURRENT_LIST_FILE}" PATH)
  include("${_dir}/dune-distrib-test-targets.cmake")
endif(dune-distrib-test_LIBRARIES)
endif(NOT dune-distrib-test_FOUND)
