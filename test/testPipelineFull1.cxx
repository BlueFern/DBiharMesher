#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectReader.h>
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
#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/centreline.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	vtkSmartPointer<vtkRescaleUnits> rescaleUnits = vtkSmartPointer<vtkRescaleUnits>::New();
	rescaleUnits->SetInputData(vesselCentreline);
	rescaleUnits->SetScale(1000); // mm to µm
	rescaleUnits->Update();

	vtkSmartPointer<vtkCentrelineResampler> centrelineSegmentSource = vtkSmartPointer<vtkCentrelineResampler>::New();
	centrelineSegmentSource->SetEdgeLength(vtkDbiharStatic::EC_AXIAL * 4);
	centrelineSegmentSource->SetInputData(rescaleUnits->GetOutput());
	centrelineSegmentSource->Update();

	vtkPolyData *resampledVesselCentreline = centrelineSegmentSource->GetOutput();

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(resampledVesselCentreline);
	scalarRadiiToVectorsFilter->Update();

	vtkSmartPointer<vtkCentrelinePartitioner> centrelinePartitioner = vtkSmartPointer<vtkCentrelinePartitioner>::New();
	centrelinePartitioner->SetInputData(scalarRadiiToVectorsFilter->GetOutput());



	vtkSmartPointer<vtkIdList> EndPoints = vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkIdList> bifurcations = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> endPointIds = vtkSmartPointer<vtkIdList>::New();

	//EndPoints->InsertNextId(180);
	//centrelinePartitioner->SetEndPoints(EndPoints);
	centrelinePartitioner->SetPartitionLength(50);
	centrelinePartitioner->Update();

	vtkPolyData *partitionedCentreline = centrelinePartitioner->GetOutput();

	vtkSmartPointer<vtkCellArray> vertexArray = vtkSmartPointer<vtkCellArray>::New();
	vtkIdType cellId = 1;
	vertexArray =  partitionedCentreline->GetVerts();
	vertexArray->GetNextCell(endPointIds);
	if (vertexArray->GetNumberOfCells() == 2)
	{
		vertexArray->GetNextCell(bifurcations);
		cellId = 2;
	}

	vtkSmartPointer<vtkAppendPolyData> fullMeshJoiner = vtkSmartPointer<vtkAppendPolyData>::New();

	bool straightSection = true;

	vtkIdType bifurcationId = -1;

	while (true)
	{
		if (cellId >= partitionedCentreline->GetNumberOfCells())
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
			vtkSmartPointer<vtkCentrelineToDbiharPatch> dbiharPatchFilter = vtkSmartPointer<vtkCentrelineToDbiharPatch>::New();
			dbiharPatchFilter->SetInputData(partitionedCentreline);

			vtkSmartPointer<vtkIdList> cellPoints = vtkSmartPointer<vtkIdList>::New();
			vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

			partitionedCentreline->GetCell(cellId, cell);
			cellPoints = cell->GetPointIds();
			for (vtkIdType id = 0; id < bifurcations->GetNumberOfIds(); id++)
			{
				vtkIdType ptId = bifurcations->GetId(id);
				vtkIdType ptPos = vtkDbiharStatic::GetPosition(cellPoints, ptId);
				if (ptPos != -1)
				{
					bifurcationId = ptPos;
					break;
				}
			}

			dbiharPatchFilter->SetNumberOfRadialQuads(28);
			dbiharPatchFilter->SetSpineId(cellId);

			dbiharPatchFilter->SetBifurcationId(bifurcationId);

			//dbiharPatchFilter->SetShowProgress(true);
			dbiharPatchFilter->SetArchDerivScale(3.2);
			dbiharPatchFilter->SetEdgeDerivScale(4.0);
			dbiharPatchFilter->Update();

			lengths[i] = partitionedCentreline->GetCell(cellId)->GetNumberOfPoints();
			appendPoints->AddInputData(dbiharPatchFilter->GetOutput());

			cellId++;

			if (i == 1) //check if just straight, if so break.
			{
				vtkSmartPointer<vtkGenericCell> cell1 = vtkSmartPointer<vtkGenericCell>::New();
				vtkSmartPointer<vtkGenericCell> cell2 = vtkSmartPointer<vtkGenericCell>::New();
				partitionedCentreline->GetCell(cellId - 2, cell1);

				previousPoints = cell1->GetPoints();
				partitionedCentreline->GetCell(cellId - 1, cell2);
				points = cell2->GetPoints();
				int tmp = points->GetNumberOfPoints() - 1;
				int tmp2 = previousPoints->GetNumberOfPoints() - 1;
				previousPoints->GetPoint(0, p0);
				points->GetPoint(tmp, p1);

				if (p0[0] == p1[0] && p0[1] == p1[1] && p0[2] == p1[2])
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
		fullMeshJoiner->AddInputData(pointsToMeshFilter->GetOutput());

	}
	fullMeshJoiner->Update();

	vtkDbiharStatic::ShowPolyData(fullMeshJoiner->GetOutput());

	vtkDbiharStatic::WritePolyData(fullMeshJoiner->GetOutput(), "quadMeshFull.vtp");

#if 0 // Very Expensive to run.

	vtkSmartPointer<vtkSubdivideMeshDynamic> subdivideMeshDynamic = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	subdivideMeshDynamic->SetInputData(appendPolyData->GetOutput());
	subdivideMeshDynamic->SetHeight(vtkDbiharStatic::EC_CIRC);
	subdivideMeshDynamic->SetLength(vtkDbiharStatic::EC_AXIAL);
	subdivideMeshDynamic->Update();

	writePolyData(subdivideMeshDynamic->GetOutput(), "ECquadMeshFull.vtp");

	vtkSmartPointer<vtkSubdivideMeshDynamic> subdivideMeshDynamic2 = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	subdivideMeshDynamic2->SetInputData(appendPolyData->GetOutput());
	subdivideMeshDynamic2->SetHeight(vtkDbiharStatic::SMC_CIRC);
	subdivideMeshDynamic2->SetLength(vtkDbiharStatic::SMC_AXIAL);
	subdivideMeshDynamic2->Update();

	writePolyData(subdivideMeshDynamic2->GetOutput(), "SMCquadMeshFull.vtp");
#endif

	vtkSmartPointer<vtkAppendPolyData> capJoiner = vtkSmartPointer<vtkAppendPolyData>::New();
	capJoiner->AddInputData(fullMeshJoiner->GetOutput());

	for (int i = 0; i < endPointIds->GetNumberOfIds(); i++)
	{
		vtkSmartPointer<vtkSkipSegmentFilter> skipSegmentFilter = vtkSmartPointer<vtkSkipSegmentFilter>::New();
		skipSegmentFilter->SetInputData(scalarRadiiToVectorsFilter->GetOutput());

		// Process only inlet on first iteration.
		skipSegmentFilter->SetInlet(i == 0);
		// Process only outlet on consecutive iterations.

		skipSegmentFilter->SetOutlet(i != 0);
		skipSegmentFilter->SetSkipSize(10);
		skipSegmentFilter->SetPointId(endPointIds->GetId(i));
		skipSegmentFilter->SetNumberOfRadialQuads(28);
		skipSegmentFilter->Update();
		capJoiner->AddInputData(skipSegmentFilter->GetOutput());

		vtkSmartPointer<vtkEndCapFilter> endCapFilter = vtkSmartPointer<vtkEndCapFilter>::New();
		endCapFilter->SetInputData(skipSegmentFilter->GetOutput());
		endCapFilter->Update();
		capJoiner->AddInputData(endCapFilter->GetOutput());
	}
	capJoiner->Update();

	vtkSmartPointer<vtkAppendPolyData> triMeshJoiner = vtkSmartPointer<vtkAppendPolyData>::New();
	triMeshJoiner->AddInputData(capJoiner->GetOutput());

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(fullMeshJoiner->GetOutput());
	triangleFilter->Update();
	triMeshJoiner->AddInputData(triangleFilter->GetOutput());

	triMeshJoiner->Update();

	vtkDbiharStatic::WriteStlData(triMeshJoiner->GetOutput(), "triMeshWithCapsFull.stl");

	return EXIT_SUCCESS;
}
