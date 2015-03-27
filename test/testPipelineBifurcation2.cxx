#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSTLWriter.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkCell.h>
#include <vtkGenericCell.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkAppendPoints.h>
#include <vtkAppendPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkCellArray.h>
#include <vtkGenericCell.h>

#include "wrapDbiharConfig.h"
#include "vtkDbiharStatic.h"
#include "vtkRescaleUnits.h"
#include "vtkCentrelinePartitioner.h"
#include "vtkScalarRadiiToVectorsFilter.h"
#include "vtkCentrelineToDbiharPatch.h"
#include "vtkPointsToMeshFilter.h"
#include "vtkSubdivideMeshDynamic.h"
#include "vtkSkipSegmentFilter.h"
#include "vtkEndCapFilter.h"
#include "vtkCentrelineResampler.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> centrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	centrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/centreline.vtk").c_str());
	centrelineReader->Update();

	vtkPolyData *centreline = vtkPolyData::SafeDownCast(centrelineReader->GetOutput());

	vtkSmartPointer<vtkRescaleUnits> rescaleUnits = vtkSmartPointer<vtkRescaleUnits>::New();
	rescaleUnits->SetInputData(centreline);
	rescaleUnits->SetScale(1000); // Rescale from mm to µm.
	rescaleUnits->Update();

	vtkSmartPointer<vtkCentrelineResampler> centrelineSegmentSource = vtkSmartPointer<vtkCentrelineResampler>::New();
	centrelineSegmentSource->SetEdgeLength(vtkDbiharStatic::EC_AXIAL * 4);
	centrelineSegmentSource->SetInputData(rescaleUnits->GetOutput());
	centrelineSegmentSource->Update();

	vtkPolyData *resampledCentreline = centrelineSegmentSource->GetOutput();

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(resampledCentreline);
	scalarRadiiToVectorsFilter->Update();

	vtkSmartPointer<vtkIdList> endPoints = vtkSmartPointer<vtkIdList>::New();
	endPoints->InsertNextId(124);
	endPoints->InsertNextId(184);
	endPoints->InsertNextId(426);

	vtkSmartPointer<vtkCentrelinePartitioner> centrelinePartitioner = vtkSmartPointer<vtkCentrelinePartitioner>::New();
	centrelinePartitioner->SetInputData(scalarRadiiToVectorsFilter->GetOutput());
	centrelinePartitioner->SetEndPoints(endPoints);
	centrelinePartitioner->SetPartitionLength(50);
	centrelinePartitioner->Update();

	vtkPolyData *partitionedCentreline = centrelinePartitioner->GetOutput();

	vtkSmartPointer<vtkIdList> endPointIds = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> bifurcations = vtkSmartPointer<vtkIdList>::New();

	// First cell (or two) are polyVertex cells describing endpoints (and bifurcations if applicable).
	int offset = 1;

	vtkSmartPointer<vtkCellArray> vertexArray = vtkSmartPointer<vtkCellArray>::New();
	vertexArray =  partitionedCentreline->GetVerts();
	vertexArray->GetNextCell(endPointIds);
	if (vertexArray->GetNumberOfCells() == 2)
	{
		vertexArray->GetNextCell(bifurcations);
		offset = 2;
	}




#if 0
	writePolyData(partitionedCentreline, "tmpCentreline2.vtp");
	std::cout << "Exiting " << __FUNCTION__ << " in " << __FILE__ << ":" << __LINE__ << std::endl;
	exit(EXIT_FAILURE);
#endif

	double lengths[3] = {0.0};
	vtkSmartPointer<vtkAppendPoints> appendPoints = vtkSmartPointer<vtkAppendPoints>::New();
	int base = 0 + offset;

	vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

	partitionedCentreline->GetCell(base, cell);
	cellPoints = cell->GetPointIds();

	vtkIdType bifurcationId = cellPoints->IsId(bifurcations->GetId(0));

	// Working with centreline partitions 3, 4, 5.
	for (int i = base; i < base + 3; i++)
	{
		vtkSmartPointer<vtkCentrelineToDbiharPatch> dbiharPatchFilter = vtkSmartPointer<vtkCentrelineToDbiharPatch>::New();
		dbiharPatchFilter->SetInputData(partitionedCentreline);
		dbiharPatchFilter->SetNumberOfRadialQuads(28);
		dbiharPatchFilter->SetSpineId(i);
		dbiharPatchFilter->SetBifurcationId(bifurcationId);
		dbiharPatchFilter->SetArchDerivScale(3.2);
		dbiharPatchFilter->SetEdgeDerivScale(4.0);
		dbiharPatchFilter->Update();

		lengths[i - base] = partitionedCentreline->GetCell(i)->GetNumberOfPoints();
		appendPoints->AddInputData(dbiharPatchFilter->GetOutput());
	}

	appendPoints->Update();

	vtkSmartPointer<vtkUnsignedIntArray> dimensions = vtkSmartPointer<vtkUnsignedIntArray>::New();
	dimensions->InsertNextValue(28);

	// Solving simultaneous equations for each branch sections length.
	dimensions->InsertNextValue((lengths[0] - lengths[1] + lengths[2]) / 2);
	dimensions->InsertNextValue((lengths[0] + lengths[1] - lengths[2]) / 2);
	dimensions->InsertNextValue((-lengths[0] + lengths[1] + lengths[2]) / 2);

	vtkSmartPointer<vtkPointsToMeshFilter> pointsToMeshFilter = vtkSmartPointer<vtkPointsToMeshFilter>::New();
	pointsToMeshFilter->SetInputData(appendPoints->GetOutput());
	pointsToMeshFilter->SetDimensions(dimensions);
	pointsToMeshFilter->Update();

	vtkDbiharStatic::ShowPolyData(pointsToMeshFilter->GetOutput());

	vtkSmartPointer<vtkXMLPolyDataWriter> quadMeshWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	quadMeshWriter->SetInputData(pointsToMeshFilter->GetOutput());
	quadMeshWriter->SetFileName("quadMeshBifurcation2.vtp");
	quadMeshWriter->Write();

#if 0
	// Very expensive to run.
	vtkSmartPointer<vtkSubdivideMeshDynamic> dynamicECMesher = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	dynamicECMesher->SetInputData(pointsToMeshFilter->GetOutput());
	dynamicECMesher->SetHeight(vtkDbiharStatic::EC_CIRC);
	dynamicECMesher->SetLength(vtkDbiharStatic::EC_AXIAL);
	dynamicECMesher->Update();

	vtkSmartPointer<vtkXMLPolyDataWriter> ECMeshWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	ECMeshWriter->SetInputData(dynamicECMesher->GetOutput());
	ECMeshWriter->SetFileName("ECquadMeshBifurcation2.vtp");
	ECMeshWriter->Write();

	vtkSmartPointer<vtkSubdivideMeshDynamic> dynamicSMCMesher = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	dynamicSMCMesher->SetInputData(pointsToMeshFilter->GetOutput());
	dynamicSMCMesher->SetHeight(vtkDbiharStatic::SMC_CIRC);
	dynamicSMCMesher->SetLength(vtkDbiharStatic::SMC_AXIAL);
	dynamicSMCMesher->Update();

	vtkSmartPointer<vtkXMLPolyDataWriter> SMCMeshWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	SMCMeshWriter->SetInputData(dynamicSMCMesher->GetOutput());
	SMCMeshWriter->SetFileName("SMCquadMeshBifurcation2.vtp");
	SMCMeshWriter->Write();
#endif

	// Create and assemble the triangulated mesh, including the triangulated quad mesh
	// along with skip segments and end caps. This is just for visualisation.
	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(pointsToMeshFilter->GetOutput());
	triangleFilter->Update();

	vtkSmartPointer<vtkAppendPolyData> appendTriMesh = vtkSmartPointer<vtkAppendPolyData>::New();
	appendTriMesh->AddInputData(triangleFilter->GetOutput());

	// Working with centreline partitions 3, 4, 5.
	for (int i = 0; i < endPointIds->GetNumberOfIds(); i++)
	{
		vtkSmartPointer<vtkSkipSegmentFilter> skipSegmentFilter = vtkSmartPointer<vtkSkipSegmentFilter>::New();
		skipSegmentFilter->SetInputData(scalarRadiiToVectorsFilter->GetOutput());

		// Only inlet on first iteration.
		skipSegmentFilter->SetInlet(i == 0);
		skipSegmentFilter->SetOutlet(i != 0);

		skipSegmentFilter->SetSkipSize(10);
		skipSegmentFilter->SetPointId(endPointIds->GetId(i));
		skipSegmentFilter->SetNumberOfRadialQuads(28);
		skipSegmentFilter->Update();

		appendTriMesh->AddInputData(skipSegmentFilter->GetOutput());

		vtkSmartPointer<vtkEndCapFilter> endCapFilter = vtkSmartPointer<vtkEndCapFilter>::New();
		endCapFilter->SetInputData(skipSegmentFilter->GetOutput());
		endCapFilter->Update();

		appendTriMesh->AddInputData(endCapFilter->GetOutput());
	}
	appendTriMesh->Update();

	vtkSmartPointer<vtkSTLWriter> triMeshWriter = vtkSmartPointer<vtkSTLWriter>::New();
	triMeshWriter->SetInputData(appendTriMesh->GetOutput());
	triMeshWriter->SetFileName("triMeshWithCapsBifurcation2.stl");
	triMeshWriter->Write();

	return EXIT_SUCCESS;
}
