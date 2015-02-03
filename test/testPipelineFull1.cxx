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
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	int check1 = resampledVesselCentreline->GetNumberOfCells();
	int lastId = 0;
	int totalNumberOfPoints = 0;

	endPointIds->InsertNextId(0);
	for (int i = 0; i < resampledVesselCentreline->GetNumberOfCells(); i++)
	{
		vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
		resampledVesselCentreline->GetCell(i, cell);
		int localEndId = cell->GetPointId(cell->GetNumberOfPoints() - 1);

		resampledVesselCentreline->GetPointCells(localEndId, startingCell);
		if (startingCell->GetNumberOfIds() == 1)
		{
			endPointIds->InsertNextId(localEndId);
		}
	}

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(resampledVesselCentreline);
	scalarRadiiToVectorsFilter->Update();

	vtkSmartPointer<vtkCentrelinePartitioner> centrelinePartitioner = vtkSmartPointer<vtkCentrelinePartitioner>::New();
	centrelinePartitioner->SetInputData(scalarRadiiToVectorsFilter->GetOutput());
	centrelinePartitioner->SetPartitionLength(100);
	centrelinePartitioner->Update();

	vtkPolyData *partitionedCentreline = centrelinePartitioner->GetOutput();

	vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();

	bool straightSection = true;
	int k = 0;

	while (true)
	{
		if (k >= partitionedCentreline->GetNumberOfCells())
		{
			break;
		}

		vtkSmartPointer<vtkAppendPoints> appendPoints = vtkSmartPointer<vtkAppendPoints>::New();
		vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPoints> previousPoints = vtkSmartPointer<vtkPoints>::New();
		double lengths[3], p0[3], p1[3] = {0.0};
		straightSection = false;

		for (int i = 0; i < 3; i++)
		{
			vtkSmartPointer<vtkCentrelineToDbiharPatch> test = vtkSmartPointer<vtkCentrelineToDbiharPatch>::New();
			test->SetInputData(partitionedCentreline);
			test->SetNumberOfRadialQuads(28);
			test->SetSpineId(k);
			test->Update();

			lengths[i] = partitionedCentreline->GetCell(k)->GetNumberOfPoints();
			appendPoints->AddInputData(test->GetOutput());

			k++;

			if (i == 1) //check if just straight, if so break.
			{
				vtkSmartPointer<vtkGenericCell> cell1 = vtkSmartPointer<vtkGenericCell>::New();
				vtkSmartPointer<vtkGenericCell> cell2 = vtkSmartPointer<vtkGenericCell>::New();
				partitionedCentreline->GetCell(k - 2, cell1);

				previousPoints = cell1->GetPoints();
				partitionedCentreline->GetCell(k - 1, cell2);
				points = cell2->GetPoints();
				int tmp = points->GetNumberOfPoints() - 1;
				int tmp2 = previousPoints->GetNumberOfPoints() - 1;
				previousPoints->GetPoint(0, p0);
				points->GetPoint(tmp, p1);

				if (p0[0] == p1[0] && p0[1] == p1[1] && p0[2] == p1[2]) // TODO: Check
				{
					straightSection = true;
					break;

				}
			}
		}

		appendPoints->Update();
		vtkSmartPointer<vtkUnsignedIntArray> dimensions = vtkSmartPointer<vtkUnsignedIntArray>::New();
		dimensions->InsertNextValue(28);

		if (straightSection)
		{
			dimensions->InsertNextValue(lengths[0] - 1);
		}
		else
		{
			// Solving simultaneous equations for each branch sections length
			dimensions->InsertNextValue((lengths[0] - lengths[1] + lengths[2]) / 2);
			dimensions->InsertNextValue((lengths[0] + lengths[1] - lengths[2]) / 2);
			dimensions->InsertNextValue((-lengths[0] + lengths[1] + lengths[2]) / 2);
		}

		vtkSmartPointer<vtkPointsToMeshFilter> pointsToMeshFilter = vtkSmartPointer<vtkPointsToMeshFilter>::New();
		pointsToMeshFilter->SetInputData(appendPoints->GetOutput());
		pointsToMeshFilter->SetDimensions(dimensions);
		pointsToMeshFilter->Update();
		appendPolyData->AddInputData(pointsToMeshFilter->GetOutput());

	}
	appendPolyData->Update();



	vtkSmartPointer<vtkXMLPolyDataWriter> writer0 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer0->SetInputData(appendPolyData->GetOutput());
	writer0->SetFileName("quadMeshFull.vtp");
	writer0->Write();

#if 0 // Very Expensive to run.

	vtkSmartPointer<vtkSubdivideMeshDynamic> subdivideMeshDynamic = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	subdivideMeshDynamic->SetInputData(appendPolyData->GetOutput());
	subdivideMeshDynamic->SetHeight(vtkDbiharStatic::EC_CIRC);
	subdivideMeshDynamic->SetLength(vtkDbiharStatic::EC_AXIAL);
	subdivideMeshDynamic->Update();

	vtkSmartPointer<vtkXMLPolyDataWriter> writer1 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer1->SetInputData(subdivideMeshDynamic->GetOutput());
	writer1->SetFileName("ECquadMeshFull.vtp");
	writer1->Write();


	vtkSmartPointer<vtkSubdivideMeshDynamic> subdivideMeshDynamic2 = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	subdivideMeshDynamic2->SetInputData(appendPolyData->GetOutput());
	subdivideMeshDynamic2->SetHeight(vtkDbiharStatic::SMC_CIRC);
	subdivideMeshDynamic2->SetLength(vtkDbiharStatic::SMC_AXIAL);
	subdivideMeshDynamic2->Update();

	vtkSmartPointer<vtkXMLPolyDataWriter> writer2 = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer2->SetInputData(subdivideMeshDynamic2->GetOutput());
	writer2->SetFileName("SMCquadMeshFull.vtp");
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
	writer3->SetFileName("triMeshWithCapsFull.stl");
	writer3->Write();

	return EXIT_SUCCESS;
}
