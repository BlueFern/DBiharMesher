cmake_minimum_required(VERSION 2.6)

project(WrapDbihar)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake ${CMAKE_MODULE_PATH})

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

add_library(vtkCentrelineResampler src/vtkCentrelineResampler.cxx)
target_link_libraries(vtkCentrelineResampler
	vtkDbiharStatic
	vtkFiltersHybrid)

add_library(vtkCentrelinePartitioner src/vtkCentrelinePartitioner.cxx)
target_link_libraries(vtkCentrelinePartitioner
	vtkDbiharStatic
	vtkFiltersExtraction
	vtkFiltersHybrid)

add_library(vtkDbiharPatchFilter src/vtkDbiharPatchFilter.cxx)
target_link_libraries(vtkDbiharPatchFilter
	dbihar
	vtkFiltersHybrid)

add_library(vtkScalarRadiiToVectorsFilter src/vtkScalarRadiiToVectorsFilter.cxx)
target_link_libraries(vtkScalarRadiiToVectorsFilter
	vtkDbiharStatic
	vtkFiltersHybrid)

add_library(vtkPointsToMeshFilter src/vtkPointsToMeshFilter.cxx)
target_link_libraries(vtkPointsToMeshFilter
	vtkDbiharStatic
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
	
add_library(vtkSubdivideMesh src/vtkSubdivideMesh.cxx)
target_link_libraries(vtkSubdivideMesh
	vtkFiltersHybrid)
	
add_library(vtkRescaleUnits src/vtkRescaleUnits.cxx)
target_link_libraries(vtkRescaleUnits
	vtkFiltersHybrid)
	
add_library(vtkDbiharPatchSmooth src/vtkDbiharPatchSmooth.cxx)
target_link_libraries(vtkDbiharPatchSmooth
	vtkFiltersHybrid)
	
include_directories(${PROJECT_SOURCE_DIR}/src)
include_directories(${PROJECT_SOURCE_DIR}/test)

add_executable(testWrapDbihar0
	test/testWrapDbihar0.cxx)
target_link_libraries(testWrapDbihar0
	vtkDbiharStatic
	vtkFiltersGeometry
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar1
	test/testWrapDbihar1.cxx)
target_link_libraries(testWrapDbihar1
	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar2
	test/testWrapDbihar2.cxx)
target_link_libraries(testWrapDbihar2
	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar3
	test/testWrapDbihar3.cxx)
target_link_libraries(testWrapDbihar3
	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar4
	test/testWrapDbihar4.cxx)
target_link_libraries(testWrapDbihar4
	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar5
	test/testWrapDbihar5.cxx)
target_link_libraries(testWrapDbihar5
	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar6
	test/testWrapDbihar6.cxx)
target_link_libraries(testWrapDbihar6
	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar7
	test/testWrapDbihar7.cxx)
target_link_libraries(testWrapDbihar7
    	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})
	
add_executable(testWrapDbihar8
	test/testWrapDbihar8.cxx)
target_link_libraries(testWrapDbihar8
    	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})
	
add_executable(testWrapDbihar9
	test/testWrapDbihar9.cxx)
target_link_libraries(testWrapDbihar9
    	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar10
	test/testWrapDbihar10.cxx)
target_link_libraries(testWrapDbihar10
    	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar11
	test/testWrapDbihar11.cxx)
target_link_libraries(testWrapDbihar11
    	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar12
	test/testWrapDbihar12.cxx)
target_link_libraries(testWrapDbihar12
    	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar13
	test/testWrapDbihar13.cxx)
target_link_libraries(testWrapDbihar13
    	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testWrapDbihar14
	test/testWrapDbihar14.cxx)
target_link_libraries(testWrapDbihar14
    	vtkDbiharStatic
	vtkDbiharPatchFilter
	${VTK_LIBRARIES})

add_executable(testCentrelineResampler
	src/vtkCentrelineResampler.cxx
	test/testCentrelineResampler.cxx)
target_link_libraries(testCentrelineResampler
	vtkDbiharStatic
	vtkCentrelineResampler
	${VTK_LIBRARIES})

add_executable(testScalarRadiiToVectorsFilter1
	src/vtkCentrelineResampler.cxx
	test/testScalarRadiiToVectorsFilter1.cxx)
target_link_libraries(testScalarRadiiToVectorsFilter1
	vtkDbiharStatic
	vtkRescaleUnits
	vtkScalarRadiiToVectorsFilter
	${VTK_LIBRARIES})

add_executable(testPointsToMeshFilter1
	test/testPointsToMeshFilter1.cxx)
target_link_libraries(testPointsToMeshFilter1
	vtkDbiharStatic
	vtkPointsToMeshFilter
	${VTK_LIBRARIES})

add_executable(testCentrelinePartitioner1
	src/vtkCentrelineResampler.cxx
	src/vtkCentrelinePartitioner.cxx
	test/testCentrelinePartitioner1.cxx)
target_link_libraries(testCentrelinePartitioner1
	vtkDbiharStatic
	vtkCentrelinePartitioner
	${VTK_LIBRARIES})

add_executable(testEndCapFilter1
	test/testEndCapFilter1.cxx)
target_link_libraries(testEndCapFilter1
	vtkDbiharStatic
	vtkEndCapFilter
	${VTK_LIBRARIES})

add_executable(testSkipSegmentFilter1
	src/vtkSkipSegmentFilter.cxx
	test/testSkipSegmentFilter1.cxx)
target_link_libraries(testSkipSegmentFilter1
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSkipSegmentFilter
	${VTK_LIBRARIES})
	
add_executable(testCentrelineToDbiharPatch
	src/vtkCentrelineToDbiharPatch.cxx
	test/testCentrelineToDbiharPatch.cxx)
target_link_libraries(testCentrelineToDbiharPatch
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkCentrelineToDbiharPatch
	vtkFiltersGeometry
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	${VTK_LIBRARIES})
	
add_executable(testSubdivideQuadFilter
	src/vtkSubdivideQuadFilter.cxx
	test/testSubdivideQuadFilter.cxx)
target_link_libraries(testSubdivideQuadFilter
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideQuadFilter
	${VTK_LIBRARIES})
	
add_executable(testSubdivideMeshDynamic
	src/vtkSubdivideQuadFilter.cxx
	src/vtkSubdivideMeshDynamic.cxx
	test/testSubdivideMeshDynamic.cxx)
target_link_libraries(testSubdivideMeshDynamic
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideQuadFilter
	vtkSubdivideMeshDynamic
	${VTK_LIBRARIES})
	
add_executable(testSubdivideMesh
	src/vtkSubdivideQuadFilter.cxx
	src/vtkSubdivideMesh.cxx
	test/testSubdivideMesh.cxx)
target_link_libraries(testSubdivideMesh
	vtkDbiharStatic
	vtkSubdivideQuadFilter
	vtkSubdivideMesh
	vtkFiltersGeometry
	${VTK_LIBRARIES})
	
add_executable(testRescaleUnits
	src/vtkRescaleUnits.cxx
	test/testRescaleUnits.cxx)
target_link_libraries(testRescaleUnits
	vtkDbiharStatic
	vtkRescaleUnits
	${VTK_LIBRARIES})
	
add_executable(testDbiharPatchSmooth
	src/vtkDbiharPatchSmooth.cxx
	test/testDbiharPatchSmooth.cxx)
target_link_libraries(testDbiharPatchSmooth
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkDbiharPatchFilter
	vtkRescaleUnits
	vtkDbiharPatchSmooth
	vtksys
	${VTK_LIBRARIES})
	
add_executable(testPipelineBifurcation1
	test/testPipelineBifurcation1.cxx)
target_link_libraries(testPipelineBifurcation1
	vtkDbiharStatic
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
	${VTK_LIBRARIES})

add_executable(testPipelineBifurcation2
	test/testPipelineBifurcation2.cxx)
target_link_libraries(testPipelineBifurcation2
	vtkDbiharStatic
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
	${VTK_LIBRARIES})

add_executable(testPipelineFull1
	test/testPipelineFull1.cxx)
target_link_libraries(testPipelineFull1
	vtkDbiharStatic
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
	${VTK_LIBRARIES})

add_executable(meshGen32
	src/meshGen32.cxx)
target_link_libraries(meshGen32
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})

add_executable(meshGen120
	src/meshGen120.cxx)
target_link_libraries(meshGen120
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})

add_executable(meshGen216
	src/meshGen216.cxx)
target_link_libraries(meshGen216
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})
	
add_executable(meshGen512
	src/meshGen512.cxx)
target_link_libraries(meshGen512
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})

add_executable(meshGen960
	src/meshGen960.cxx)
target_link_libraries(meshGen960
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})

add_executable(meshGen2000
	src/meshGen2000.cxx)
target_link_libraries(meshGen2000
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})

add_executable(meshGen4000
	src/meshGen4000.cxx)
target_link_libraries(meshGen4000
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})

add_executable(testPipelineFullc4032
	test/testPipelineFullc4032.cxx)
target_link_libraries(testPipelineFullc4032
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})

add_executable(meshGen4080
	src/meshGen4080.cxx)
target_link_libraries(meshGen4080
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	vtksys
	${VTK_LIBRARIES})
	
add_executable(meshGenSmooth4080
	src/meshGenSmooth4080.cxx)
target_link_libraries(meshGenSmooth4080
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	vtksys
	vtkDbiharPatchSmooth
	${VTK_LIBRARIES})

add_executable(meshGen8000
	src/meshGen8000.cxx)
target_link_libraries(meshGen8000
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})

add_executable(meshGen8064
	src/meshGen8064.cxx)
target_link_libraries(meshGen8064
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})
	
add_executable(meshGenSmooth8064
	src/meshGenSmooth8064.cxx)
target_link_libraries(meshGenSmooth8064
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	vtksys
	vtkDbiharPatchSmooth
	${VTK_LIBRARIES})

add_executable(meshGen16128
	src/meshGen16128.cxx)
target_link_libraries(meshGen16128
	vtkDbiharStatic
	vtkScalarRadiiToVectorsFilter
	vtkFiltersGeometry
	vtkSubdivideMesh
	vtkSubdivideQuadFilter
	vtkCentrelineToDbiharPatch
	vtkDbiharPatchFilter
	vtkCentrelinePartitioner
	vtkPointsToMeshFilter
	vtkSkipSegmentFilter
	vtkEndCapFilter
	vtkRescaleUnits
	vtkCentrelineResampler
	${VTK_LIBRARIES})

#add_executable(DumpMeshToLegacyFormat
#	util/DumpMeshToLegacyFormat.cxx)
#target_link_libraries(DumpMeshToLegacyFormat
#	vtkIOXML vtkIOLegacy)

#add_executable(scratch
#	test/showPolyData.cxx
#	test/scratch.cxx)
#target_link_libraries(scratch
#	${VTK_LIBRARIES}
#	vtkIOLegacy)
