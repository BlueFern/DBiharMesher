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

add_subdirectory(libFortran/dbihar)

foreach(vtkLib ${VTK_LIBRARIES})
	message(STATUS ${vtkLib})
endforeach(vtkLib)

add_library(vtkDbiharPatchFilter src/vtkDbiharPatchFilter.cxx)
target_link_libraries(vtkDbiharPatchFilter
	dbihar
	vtkFiltersHybrid)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/test)

message(STATUS ${VTK_LIBRARIES})

set(testLibs vtkDbiharPatchFilter vtkRenderingAnnotation vtkRenderingFreeTypeOpenGL vtkInteractionStyle vtkFiltersSources)

add_executable(testWrapDbihar1
	test/showPolyData.cxx
	test/testWrapDbihar1.cxx)
target_link_libraries(testWrapDbihar1
	${testLibs})

add_executable(testWrapDbihar2
	test/showPolyData.cxx
	test/testWrapDbihar2.cxx)
target_link_libraries(testWrapDbihar2
	${testLibs})

add_executable(testWrapDbihar3
	test/showPolyData.cxx
	test/testWrapDbihar3.cxx)
target_link_libraries(testWrapDbihar3
	${testLibs})

add_executable(testWrapDbihar4
	test/showPolyData.cxx
	test/testWrapDbihar4.cxx)
target_link_libraries(testWrapDbihar4
	${testLibs})

add_executable(testWrapDbihar5
	test/showPolyData.cxx
	test/testWrapDbihar5.cxx)
target_link_libraries(testWrapDbihar5
	${testLibs})
	