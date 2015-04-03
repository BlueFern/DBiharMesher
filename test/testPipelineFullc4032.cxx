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
#include "vtkSubdivideMesh.h"
#include "vtkSkipSegmentFilter.h"
#include "vtkEndCapFilter.h"
#include "vtkCentrelineResampler.h"
#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/../tmpData/c4032/c4032Centreline.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	vtkSmartPointer<vtkRescaleUnits> rescaleUnits = vtkSmartPointer<vtkRescaleUnits>::New();
	rescaleUnits->SetInputData(vesselCentreline);
	rescaleUnits->SetScale(1000); // Convert mm to µm.
	rescaleUnits->Update();

	vtkSmartPointer<vtkCentrelineResampler> centrelineSegmentSource = vtkSmartPointer<vtkCentrelineResampler>::New();
	centrelineSegmentSource->SetEdgeLength(vtkDbiharStatic::EC_AXIAL * 4);
	centrelineSegmentSource->SetInputData(rescaleUnits->GetOutput());
	centrelineSegmentSource->Update();

	// vtkDbiharStatic::ShowPolyData(centrelineSegmentSource->GetOutput());

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(centrelineSegmentSource->GetOutput());
	scalarRadiiToVectorsFilter->Update();

	// vtkDbiharStatic::ShowPolyData(scalarRadiiToVectorsFilter->GetOutput());

	vtkSmartPointer<vtkCentrelinePartitioner> centrelinePartitioner = vtkSmartPointer<vtkCentrelinePartitioner>::New();
	centrelinePartitioner->SetInputData(scalarRadiiToVectorsFilter->GetOutput());
	centrelinePartitioner->SetPartitionLength(100);
	centrelinePartitioner->Update();

	vtkDbiharStatic::ShowPolyData(centrelinePartitioner->GetOutput());
	// vtkDbiharStatic::WritePolyData(centrelinePartitioner->GetOutput(), "c4032_part.vtp");

	vtkPolyData *partitionedCentreline = centrelinePartitioner->GetOutput();

	vtkSmartPointer<vtkAppendPolyData> fullMeshJoiner = vtkSmartPointer<vtkAppendPolyData>::New();

	vtkSmartPointer<vtkCellArray> allSpines = partitionedCentreline->GetLines();

	// WARNING: There is a duplication of the following code in vtkCentrelineToDbiharPatch.

	// Get end points and bifurcations.
	vtkSmartPointer<vtkCellArray> verts = partitionedCentreline->GetVerts();
	verts->InitTraversal();

	// Get end points for later.
	// WARNING: the assumption here is that end points always come first.
	vtkSmartPointer<vtkIdList> endPointIds = vtkSmartPointer<vtkIdList>::New();
	verts->GetNextCell(endPointIds);

	// Get bifurcation points.
	vtkSmartPointer<vtkIdList> bifurcationIds = vtkSmartPointer<vtkIdList>::New();
	verts->GetNextCell(bifurcationIds);

	int numRadialQuads = 48;
	vtkSmartPointer<vtkAppendPoints> appendPoints;
	vtkSmartPointer<vtkUnsignedIntArray> pointsToMeshDimensions;

	// Number of inputs for the appendPoints filter.
	unsigned int inputId = 0;

	// For straight segments this is 1, for bifurcations this will become 2.
	unsigned int maxInputId = 1;

	// Iterate over all spines.
	vtkIdType spineId = 0;
	while(true)
	{
		// Find the spine indicated by the spineId.
		allSpines->InitTraversal();
		vtkSmartPointer<vtkIdList> spineIds = vtkSmartPointer<vtkIdList>::New();
		vtkIdType searchId = 0;
		while(true)
		{
			allSpines->GetNextCell(spineIds);
			if(searchId == spineId)
			{
				break;
			}
			else
			{
				searchId++;
			}
		}
		const int spineSize = spineIds->GetNumberOfIds();

		// Determine if this is a bifurcation.
		// Remember the position if this is the case.
		// TODO: Check the spine contains only one bifurcation.
		bool bifurcation = false;
		vtkIdType bifurcationPos = -1;
		for(vtkIdType bifPos = 0; bifPos < bifurcationIds->GetNumberOfIds(); bifPos++)
		{
			bifurcationPos = spineIds->IsId(bifurcationIds->GetId(bifPos));
			if(bifurcationPos != -1)
			{
				bifurcation = true;
				maxInputId = 2;
				break;
			}
		}

		// Prepare to append points.
		if(inputId == 0)
		{
			appendPoints = vtkSmartPointer<vtkAppendPoints>::New();
			pointsToMeshDimensions = vtkSmartPointer<vtkUnsignedIntArray>::New();
			pointsToMeshDimensions->InsertNextValue(numRadialQuads);
		}

		vtkSmartPointer<vtkCentrelineToDbiharPatch> dbiharPatchFilter = vtkSmartPointer<vtkCentrelineToDbiharPatch>::New();
		dbiharPatchFilter->SetInputData(partitionedCentreline);
		dbiharPatchFilter->SetNumberOfRadialQuads(numRadialQuads);
		dbiharPatchFilter->SetSpineId(spineId);
		dbiharPatchFilter->SetArchDerivScale(3.2);
		dbiharPatchFilter->SetEdgeDerivScale(4.0);
		dbiharPatchFilter->Update();

		// Set dimensions for a bifurcation.
		if(bifurcation)
		{
			pointsToMeshDimensions->InsertNextValue(bifurcationPos);
		}

		// Send output to appendPoints.
		appendPoints->AddInputData(dbiharPatchFilter->GetOutput());
		inputId++;
		spineId++;

		// All inputs to appendPoints are ready.
		if(inputId > maxInputId)
		{
			appendPoints->Update();

			// Set dimensions for a straight segment.
			if(!bifurcation)
			{
				pointsToMeshDimensions->InsertNextValue(spineSize - 1);
			}

			// Run pointsToMesh filter.
			vtkSmartPointer<vtkPointsToMeshFilter> pointsToMeshFilter = vtkSmartPointer<vtkPointsToMeshFilter>::New();
			pointsToMeshFilter->SetInputData(appendPoints->GetOutput());
			pointsToMeshFilter->SetDimensions(pointsToMeshDimensions);
			pointsToMeshFilter->Update();

			// Send output to the append filter.
			fullMeshJoiner->AddInputData(pointsToMeshFilter->GetOutput());

			inputId = 0;
			break;
		}
	}

	fullMeshJoiner->Update();

	vtkDbiharStatic::ShowPolyData(fullMeshJoiner->GetOutput());
	vtkDbiharStatic::WritePolyData(fullMeshJoiner->GetOutput(), "quadMeshFullc4032.vtp");

#if 1
	int numECs = 4;
	int numSMCs = 4;
	vtkSmartPointer<vtkSubdivideMesh> subdivideECMesh = vtkSmartPointer<vtkSubdivideMesh>::New();
	//subdivideECMesh->DebugOn();
	subdivideECMesh->SetInputData(fullMeshJoiner->GetOutput());
	subdivideECMesh->SetColumns(numECs);
	subdivideECMesh->SetRows((vtkDbiharStatic::SMC_CIRC / vtkDbiharStatic::EC_CIRC) * numSMCs);
	//subdivideECMesh->SetHeight(vtkDbiharStatic::EC_CIRC);
	//subdivideECMesh->SetLength(vtkDbiharStatic::EC_AXIAL);
	//subdivideECMesh->Print(std::cout);
	subdivideECMesh->Update();

	vtkDbiharStatic::WritePolyData(subdivideECMesh->GetOutput(), "quadMeshFullECc4032.vtp");

	vtkSmartPointer<vtkSubdivideMesh> subdivideSMCMesh = vtkSmartPointer<vtkSubdivideMesh>::New();
	subdivideSMCMesh->SetInputData(fullMeshJoiner->GetOutput());
	//subdivideSMCMesh->DebugOn();
	subdivideSMCMesh->SetColumns((vtkDbiharStatic::EC_AXIAL / vtkDbiharStatic::SMC_AXIAL) * numECs);
	subdivideSMCMesh->SetRows(numSMCs);
	//subdivideSMCMesh->SetHeight(vtkDbiharStatic::SMC_CIRC);
	//subdivideSMCMesh->SetLength(vtkDbiharStatic::SMC_AXIAL);
	//subdivideSMCMesh->Print(std::cout);
	subdivideSMCMesh->Update();

	vtkDbiharStatic::WritePolyData(subdivideSMCMesh->GetOutput(), "quadMeshFullSMCc4032.vtp");
#endif

#if 1
	vtkSmartPointer<vtkAppendPolyData> capJoiner = vtkSmartPointer<vtkAppendPolyData>::New();

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
		skipSegmentFilter->SetNumberOfRadialQuads(numRadialQuads);
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

	vtkDbiharStatic::WriteStlData(triMeshJoiner->GetOutput(), "triMeshWithCapsFullc4032.stl");
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
