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

find_package(Doxygen)

find_Package(Graphviz)
if(DOXYGEN_DOT_EXECUTABLE)
	set(HAVE_DOT "YES")
	set(CALL_GRAPH "YES")
	set(CALLER_GRAPH "YES")
else(DOXYGEN_DOT_EXECUTABLE)
	set(HAVE_DOT "NO")
	set(CALL_GRAPH "NO")
	set(CALLER_GRAPH "NO")
endif(DOXYGEN_DOT_EXECUTABLE)

set(DOXY_INPUT_DIR ${PROJECT_SOURCE_DIR}/src)
set(DOXY_BASE_DIR ${PROJECT_SOURCE_DIR}/doc)
set(DOXY_OUTPUT_DIR ${PROJECT_BINARY_DIR}/doc)
file(MAKE_DIRECTORY ${DOXY_OUTPUT_DIR})

set(DOXY_API_LINK_TAGS ${DOXY_OUTPUT_DIR}/DoxyLink.tag)

configure_file(
	${PROJECT_SOURCE_DIR}/doc/Doxyfile.in
	${PROJECT_BINARY_DIR}/Doxyfile
	@ONLY)

if(DOXYGEN_FOUND)
	add_custom_target(doxyDocs
	COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile
	WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
	COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)


foreach(vtkLib ${VTK_LIBRARIES})
	message(STATUS ${vtkLib})
endforeach(vtkLib)

add_library(vtkDbiharStatic src/vtkDbiharStatic.cxx)
target_link_libraries(vtkDbiharStatic
	vtkCommonCore)

add_library(vtkDbiharPatchFilter src/vtkDbiharPatchFilter.cxx)
target_link_libraries(vtkDbiharPatchFilter
	dbihar
	vtkFiltersHybrid)

add_library(vtkScalarRadiiToVectorsFilter src/vtkScalarRadiiToVectorsFilter.cxx)
target_link_libraries(vtkScalarRadiiToVectorsFilter
	vtkDbiharStatic
	vtkFiltersHybrid)

add_library(vtkEndPointIdsToDbiharPatchFilter src/vtkEndPointIdsToDbiharPatchFilter.cxx)
target_link_libraries(vtkEndPointIdsToDbiharPatchFilter
	vtkFiltersHybrid)

add_library(vtkPointsToMeshFilter src/vtkPointsToMeshFilter.cxx)
target_link_libraries(vtkPointsToMeshFilter
	vtkDbiharStatic
	vtkFiltersHybrid)

add_library(vtkCentrelinePartitioner src/vtkCentrelinePartitioner.cxx)
target_link_libraries(vtkCentrelinePartitioner
	vtkDbiharStatic
	vtkFiltersExtraction
	vtkFiltersHybrid)

add_library(vtkEndCapFilter src/vtkEndCapFilter.cxx)
target_link_libraries(vtkEndCapFilter
	vtkFiltersHybrid
	vtkFiltersGeometry
	vtkDbiharPatchFilter)

add_library(vtkSkipSegmentFilter src/vtkSkipSegmentFilter.cxx)
target_link_libraries(vtkSkipSegmentFilter
	vtkDbiharStatic
	vtkFiltersHybrid)
	
add_library(vtkCentrelineToDbiharPatch src/vtkCentrelineToDbiharPatch.cxx)
target_link_libraries(vtkCentrelineToDbiharPatch
	vtkDbiharStatic
	vtkFiltersHybrid)
	
add_library(vtkSubdivideQuadFilter src/vtkSubdivideQuadFilter.cxx)
target_link_libraries(vtkSubdivideQuadFilter
	vtkFiltersHybrid)
	
add_library(vtkSubdivideMeshDynamic src/vtkSubdivideMeshDynamic.cxx)
target_link_libraries(vtkSubdivideMeshDynamic
	vtkFiltersHybrid)
	
add_library(vtkRescaleUnits src/vtkRescaleUnits.cxx)
target_link_libraries(vtkRescaleUnits
	vtkFiltersHybrid)
	
add_library(vtkCentrelineResampler src/vtkCentrelineResampler.cxx)
target_link_libraries(vtkCentrelineResampler
	vtkFiltersHybrid)

include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/test)

set(testLibs vtkIOXML vtkIOLegacy vtkIOImage vtkRenderingAnnotation vtkRenderingFreeTypeOpenGL vtkInteractionStyle vtkFiltersSources)

add_executable(testWrapDbihar0
	test/showPolyData.cxx
	test/testWrapDbihar0.cxx)
target_link_libraries(testWrapDbihar0
	vtkFiltersGeometry
	vtkDbiharPatchFilter
	${testLibs})

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
	vtkDbiharStatic
	vtkFiltersGeometry
	vtkScalarRadiiToVectorsFilter
	vtkEndPointIdsToDbiharPatchFilter
	vtkDbiharPatchFilter
	${testLibs})

add_executable(testEndPointIdsToDbiharPatchFilter2
	src/vtkCentrelineData.cxx
	test/showPolyData.cxx
	test/testEndPointIdsToDbiharPatchFilter2.cxx)
target_link_libraries(testEndPointIdsToDbiharPatchFilter2
	vtkDbiharStatic
	vtkFiltersGeometry
	vtkScalarRadiiToVectorsFilter
	vtkEndPointIdsToDbiharPatchFilter
	vtkDbiharPatchFilter
	${testLibs})

add_executable(testPointsToMeshFilter1
	test/testPointsToMeshFilter.cxx
	test/showPolyData.cxx)
target_link_libraries(testPointsToMeshFilter1
	vtkPointsToMeshFilter
	${testLibs})

add_executable(testCentrelinePartitioner1
	src/vtkCentrelineData.cxx
	src/vtkCentrelinePartitioner.cxx
	test/testCentrelinePartitioner1.cxx)
target_link_libraries(testCentrelinePartitioner1
	vtkCentrelinePartitioner
	${testLibs})

add_executable(testEndCapFilter1
	test/showPolyData.cxx
	test/testEndCapFilter1.cxx)
target_link_libraries(testEndCapFilter1
	vtkEndCapFilter
	${testLibs})

add_executable(testSkipSegmentFilter1
	src/vtkSkipSegmentFilter.cxx
	test/testSkipSegmentFilter1.cxx)
target_link_libraries(testSkipSegmentFilter1
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSkipSegmentFilter
	${testLibs})
	
add_executable(testCentrelineToDbiharPatch
	src/vtkCentrelineToDbiharPatch.cxx
	test/showPolyData.cxx
	test/testCentrelineToDbiharPatch.cxx)
target_link_libraries(testCentrelineToDbiharPatch
	vtkScalarRadiiToVectorsFilter
	vtkCentrelineToDbiharPatch
	vtkFiltersGeometry
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	${testLibs})
	
add_executable(testSubdivideQuadFilter
	src/vtkSubdivideQuadFilter.cxx
	test/testSubdivideQuadFilter.cxx)
target_link_libraries(testSubdivideQuadFilter
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideQuadFilter
	${testLibs})
	
add_executable(testSubdivideMeshStatic
	src/vtkSubdivideQuadFilter.cxx
	test/testSubdivideMeshStatic.cxx)
target_link_libraries(testSubdivideMeshStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideQuadFilter
	${testLibs})
	
add_executable(testSubdivideMeshDynamic
	src/vtkSubdivideQuadFilter.cxx
	src/vtkSubdivideMeshDynamic.cxx
	test/testSubdivideMeshDynamic.cxx)
target_link_libraries(testSubdivideMeshDynamic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideQuadFilter
	vtkSubdivideMeshDynamic
	${testLibs})
	
add_executable(testRescaleUnits
	src/vtkRescaleUnits.cxx
	test/testRescaleUnits.cxx)
target_link_libraries(testRescaleUnits
	vtkRescaleUnits
	${testLibs})
	
add_executable(testPipelineBifurcation1
	test/showPolyData.cxx
	test/testPipelineBifurcation1.cxx)
target_link_libraries(testPipelineBifurcation1
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMeshDynamic
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${testLibs})
	
add_executable(testPipelineFull1
	test/showPolyData.cxx
	test/testPipelineFull1.cxx)
target_link_libraries(testPipelineFull1
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMeshDynamic
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${testLibs})

#add_executable(scratch
#	test/showPolyData.cxx
#	test/scratch.cxx)
#target_link_libraries(scratch
#	${testLibs}
#	vtkIOLegacy)
