#include <map>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkMath.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>
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
#include <vtkXMLStructuredGridReader.h>
#include "showPolyData.h"
#include "vtkDbiharStatic.h"
#include "vtkRescaleUnits.h"
#include "vtkCentrelineData.h"
#include "vtkCentrelinePartitioner.h"
#include "vtkScalarRadiiToVectorsFilter.h"
#include "vtkCentrelineToDbiharPatch.h"
#include "vtkPointsToMeshFilter.h"
#include "vtkSubdivideMeshDynamic.h"
#include "vtkSkipSegmentFilter.h"
#include "vtkEndCapFilter.h"
#include "vtkCentrelineData.h"
#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/centreline.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	vtkSmartPointer<vtkRescaleUnits> rescaleUnits = vtkSmartPointer<vtkRescaleUnits>::New();
	rescaleUnits->SetInputData(vesselCentreline);
	rescaleUnits->SetScale(1000); // mm to µ
	rescaleUnits->Update();

	vtkSmartPointer<vtkCentrelineData> centrelineSegmentSource = vtkSmartPointer<vtkCentrelineData>::New();
	centrelineSegmentSource->DebugOn();
	centrelineSegmentSource->SetEdgeLength(vtkDbiharStatic::EC_AXIAL * 4);
	centrelineSegmentSource->SetCentrelineData(rescaleUnits->GetOutput());

	vtkPolyData *resampledVesselCentreline = centrelineSegmentSource->GetOutput();
	vtkSmartPointer<vtkDoubleArray> scalars = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkIdList> endPointIds = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> startingCell = vtkSmartPointer<vtkIdList>::New();

	int check1 = resampledVesselCentreline->GetNumberOfCells();
	int lastId = 0;
	int totalNumberOfPoints = 0;

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(resampledVesselCentreline);
	scalarRadiiToVectorsFilter->Update();

	vtkSmartPointer<vtkCentrelinePartitioner> centrelinePartitioner = vtkSmartPointer<vtkCentrelinePartitioner>::New();
	centrelinePartitioner->SetInputData(scalarRadiiToVectorsFilter->GetOutput());
	centrelinePartitioner->SetPartitionLength(100);
	centrelinePartitioner->Update();

	vtkPolyData *partitionedCentreline = centrelinePartitioner->GetOutput();

	vtkSmartPointer<vtkGenericCell> cell1 = vtkSmartPointer<vtkGenericCell>::New();
	vtkSmartPointer<vtkGenericCell> cell2 = vtkSmartPointer<vtkGenericCell>::New();

	partitionedCentreline->GetCell(3, cell1);

	endPointIds->InsertNextId(cell1->GetPointId(0));

	partitionedCentreline->GetCell(4, cell2);

	endPointIds->InsertNextId(cell2->GetPointId(0));
	endPointIds->InsertNextId(cell2->GetPointId(cell2->GetNumberOfPoints() - 1));


	vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();

	vtkSmartPointer<vtkAppendPoints> appendPoints = vtkSmartPointer<vtkAppendPoints>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> previousPoints = vtkSmartPointer<vtkPoints>::New();
	double lengths[3] = {0.0};

	for (int i = 3; i < 6; i++)
	{
		vtkSmartPointer<vtkCentrelineToDbiharPatch> test = vtkSmartPointer<vtkCentrelineToDbiharPatch>::New();
		test->SetInputData(partitionedCentreline);
		test->SetNumberOfRadialQuads(28);
		test->SetSpineId(i);
		test->Update();

		lengths[i] = partitionedCentreline->GetCell(i)->GetNumberOfPoints();
		appendPoints->AddInputData(test->GetOutput());
	}

	appendPoints->Update();
	vtkSmartPointer<vtkUnsignedIntArray> dimensions = vtkSmartPointer<vtkUnsignedIntArray>::New();
	dimensions->InsertNextValue(28);

	// Solving simultaneous equations for each branch sections length
	dimensions->InsertNextValue((lengths[3] - lengths[4] + lengths[5]) / 2);
	dimensions->InsertNextValue((lengths[3] + lengths[4] - lengths[5]) / 2);
	dimensions->InsertNextValue((-lengths[3] + lengths[4] + lengths[5]) / 2);

	vtkSmartPointer<vtkPointsToMeshFilter> pointsToMeshFilter = vtkSmartPointer<vtkPointsToMeshFilter>::New();
	pointsToMeshFilter->SetInputData(appendPoints->GetOutput());
	pointsToMeshFilter->SetDimensions(dimensions);
	pointsToMeshFilter->Update();
	appendPolyData->AddInputData(pointsToMeshFilter->GetOutput());

	appendPolyData->Update();

	vtkSmartPointer<vtkXMLPolyDataWriter> writer0 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer0->SetInputData(appendPolyData->GetOutput());
	writer0->SetFileName("quadMeshBifurcation.vtp");
	writer0->Write();

#if 0 // Very Expensive to run.

	vtkSmartPointer<vtkSubdivideMeshDynamic> subdivideMeshDynamic = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	subdivideMeshDynamic->SetInputData(appendPolyData->GetOutput());
	subdivideMeshDynamic->SetHeight(vtkDbiharStatic::EC_CIRC);
	subdivideMeshDynamic->SetLength(vtkDbiharStatic::EC_AXIAL);
	subdivideMeshDynamic->Update();

	vtkSmartPointer<vtkXMLPolyDataWriter> writer1 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer1->SetInputData(subdivideMeshDynamic->GetOutput());
	writer1->SetFileName("ECquadMeshBifurcation.vtp");
	writer1->Write();


	vtkSmartPointer<vtkSubdivideMeshDynamic> subdivideMeshDynamic2 = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	subdivideMeshDynamic2->SetInputData(appendPolyData->GetOutput());
	subdivideMeshDynamic2->SetHeight(vtkDbiharStatic::SMC_CIRC);
	subdivideMeshDynamic2->SetLength(vtkDbiharStatic::SMC_AXIAL);
	subdivideMeshDynamic2->Update();

	vtkSmartPointer<vtkXMLPolyDataWriter> writer2 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer2->SetInputData(subdivideMeshDynamic2->GetOutput());
	writer2->SetFileName("SMCquadMeshBifurcation.vtp");
	writer2->Write();

#endif

	vtkSmartPointer<vtkAppendPolyData> appendPolyData2 = vtkSmartPointer<vtkAppendPolyData>::New();
	vtkSmartPointer<vtkAppendPolyData> appendPolyData3 = vtkSmartPointer<vtkAppendPolyData>::New();

	appendPolyData2->AddInputData(appendPolyData->GetOutput());

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
		appendPolyData2->AddInputData(skipSegmentFilter->GetOutput());

		vtkSmartPointer<vtkEndCapFilter> endCapFilter = vtkSmartPointer<vtkEndCapFilter>::New();
		endCapFilter->SetInputData(skipSegmentFilter->GetOutput());
		endCapFilter->Update();
		appendPolyData2->AddInputData(endCapFilter->GetOutput());
	}
	appendPolyData2->Update();
	appendPolyData3->AddInputData(appendPolyData2->GetOutput());

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(appendPolyData->GetOutput());
	triangleFilter->Update();
	appendPolyData3->AddInputData(triangleFilter->GetOutput());

	appendPolyData3->Update();

	vtkSmartPointer<vtkSTLWriter> writer3 = vtkSmartPointer<vtkSTLWriter>::New();
	writer3->SetInputData(appendPolyData3->GetOutput());
	writer3->SetFileName("triMeshWithCapsBifurcation.stl");
	writer3->Write();

	return EXIT_SUCCESS;
}
