cmake_minimum_required(VERSION 2.6)

project(WrapDbihar)

set(TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/data CACHE PATH "Location of input data for the test drivers.")
configure_file(
	${PROJECT_SOURCE_DIR}/wrapDbiharConfig.h.in
	${PROJECT_BINARY_DIR}/wrapDbiharConfig.h)
include_directories(${PROJECT_BINARY_DIR})

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

add_library(vtkScalarRadiiToVectorsFilter src/vtkScalarRadiiToVectorsFilter.cxx)
target_link_libraries(vtkScalarRadiiToVectorsFilter
	vtkFiltersHybrid)

add_library(vtkEndPointIdsToDbiharPatchFilter src/vtkEndPointIdsToDbiharPatchFilter.cxx)
target_link_libraries(vtkEndPointIdsToDbiharPatchFilter
	vtkFiltersHybrid)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/test)

set(testLibs vtkIOXML vtkIOLegacy vtkIOImage vtkRenderingAnnotation vtkRenderingFreeTypeOpenGL vtkInteractionStyle vtkFiltersSources)

add_executable(testWrapDbihar1
	test/showPolyData.cxx
	test/testWrapDbihar1.cxx)
target_link_libraries(testWrapDbihar1
	vtkDbiharPatchFilter
	${testLibs})

add_executable(testWrapDbihar2
	test/showPolyData.cxx
	test/testWrapDbihar2.cxx)
target_link_libraries(testWrapDbihar2
	vtkDbiharPatchFilter
	${testLibs})

add_executable(testWrapDbihar3
	test/showPolyData.cxx
	test/testWrapDbihar3.cxx)
target_link_libraries(testWrapDbihar3
	vtkDbiharPatchFilter
	${testLibs})

add_executable(testWrapDbihar4
	test/showPolyData.cxx
	test/testWrapDbihar4.cxx)
target_link_libraries(testWrapDbihar4
	vtkDbiharPatchFilter
	${testLibs})

add_executable(testWrapDbihar5
	test/showPolyData.cxx
	test/testWrapDbihar5.cxx)
target_link_libraries(testWrapDbihar5
	vtkDbiharPatchFilter
	${testLibs})

add_executable(testWrapDbihar6
	test/showPolyData.cxx
	test/testWrapDbihar6.cxx)
target_link_libraries(testWrapDbihar6
	vtkDbiharPatchFilter
	${testLibs})

add_executable(testCentrelineData1
	src/vtkCentrelineData.cxx
	test/showPolyData.cxx
	test/testCentrelineData1.cxx)
target_link_libraries(testCentrelineData1
	vtkDbiharPatchFilter
	${testLibs})

add_executable(testScalarRadiiToVectorsFilter1
	src/vtkCentrelineData.cxx
	test/showPolyData.cxx
	test/testScalarRadiiToVectorsFilter1.cxx)
target_link_libraries(testScalarRadiiToVectorsFilter1
	vtkDbiharPatchFilter
	vtkScalarRadiiToVectorsFilter
	${testLibs})

add_executable(testEndPointIdsToDbiharPatchFilter1
	src/vtkCentrelineData.cxx
	test/showPolyData.cxx
	test/testEndPointIdsToDbiharPatchFilter1.cxx)
target_link_libraries(testEndPointIdsToDbiharPatchFilter1
	vtkScalarRadiiToVectorsFilter
	vtkEndPointIdsToDbiharPatchFilter
	vtkDbiharPatchFilter
	${testLibs})

add_executable(testEndPointIdsToDbiharPatchFilter2
	src/vtkCentrelineData.cxx
	test/showPolyData.cxx
	test/testEndPointIdsToDbiharPatchFilter2.cxx)
target_link_libraries(testEndPointIdsToDbiharPatchFilter2
	vtkScalarRadiiToVectorsFilter
	vtkEndPointIdsToDbiharPatchFilter
	vtkDbiharPatchFilter
	${testLibs})

#add_executable(scratch
#	test/showPolyData.cxx
#	test/scratch.cxx)
#target_link_libraries(scratch
#	${testLibs}
#	vtkIOLegacy)
