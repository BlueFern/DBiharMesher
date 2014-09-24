cmake_minimum_required(VERSION 2.6)

project(WrapDbihar)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
message(STATUS ${VTK_USE_FILE})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: None Debug Release."
    FORCE)
endif(NOT CMAKE_BUILD_TYPE)

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # Location of malloc.h on OS X.
  include_directories(/usr/include/malloc/)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

enable_language(Fortran)
if(CMAKE_Fortran_COMPILER)
  message(STATUS "Using Fortran compiler: " ${CMAKE_Fortran_COMPILER})
else()
  message(FATAL_ERROR "No Fortran support. Unable to continue")
endif()

add_subdirectory(libFortran/blas)
add_subdirectory(libFortran/linpack)
add_subdirectory(libFortran/dbihar)

foreach(vtkLib ${VTK_LIBRARIES})
  message(STATUS ${vtkLib})
endforeach(vtkLib)

add_library(vtkDbiharPatchFilter src/vtkDbiharPatchFilter.cxx)
target_link_libraries(vtkDbiharPatchFilter
  dbihar
  linpack
  blas
  vtkInteractionStyle
  vtkRenderingFreeTypeOpenGL
  vtkFiltersHybrid)

#   vtkRenderingCore

include_directories (${PROJECT_SOURCE_DIR}/src)
add_executable(testWrapDbihar test/testWrapDbihar.cxx)
target_link_libraries(testWrapDbihar
  vtkDbiharPatchFilter)
