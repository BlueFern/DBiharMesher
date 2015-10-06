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

#include <sstream>
#define SSTR( x ) dynamic_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()

// Run this code in the directory containing the centreline input file.

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	// Half a ring.
	int numRadialQuads = 20;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName("c8000Centreline.vtk");
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
	centrelinePartitioner->SetPartitionLength(300);
	centrelinePartitioner->Update();

	// vtkDbiharStatic::ShowPolyData(centrelinePartitioner->GetOutput());

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
		// TODO: Surely there's an easier way for this traversal.
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

		// Prepare to append points.
		if(inputId == 0)
		{
			appendPoints = vtkSmartPointer<vtkAppendPoints>::New();
			pointsToMeshDimensions = vtkSmartPointer<vtkUnsignedIntArray>::New();
			pointsToMeshDimensions->InsertNextValue(numRadialQuads);
			pointsToMeshDimensions->InsertNextValue(spineSize - 1);
		}

		vtkSmartPointer<vtkCentrelineToDbiharPatch> dbiharPatchFilter = vtkSmartPointer<vtkCentrelineToDbiharPatch>::New();
		dbiharPatchFilter->SetInputData(partitionedCentreline);
		dbiharPatchFilter->SetNumberOfRadialQuads(numRadialQuads);
		dbiharPatchFilter->SetSpineId(spineId);
		dbiharPatchFilter->SetArchDerivScale(3.2);
		dbiharPatchFilter->SetEdgeDerivScale(4.0);
		dbiharPatchFilter->Update();

		// Send output to appendPoints.
		appendPoints->AddInputData(dbiharPatchFilter->GetOutput());
		inputId++;
		spineId++;

		// All inputs to appendPoints are ready.
		if(inputId > maxInputId)
		{
			appendPoints->Update();

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
	// vtkDbiharStatic::WritePolyData(fullMeshJoiner->GetOutput(), (prefix + "quadMeshFullc8000.vtp").c_str());
	vtkDbiharStatic::WritePolyData(fullMeshJoiner->GetOutput(), "quadMeshFullc8000.vtp");

	int numECs = 4;
	int numSMCs = 4;
	vtkSmartPointer<vtkSubdivideMesh> subdivideECMesh = vtkSmartPointer<vtkSubdivideMesh>::New();
	subdivideECMesh->SetInputData(fullMeshJoiner->GetOutput());
	subdivideECMesh->SetRows(numECs);
	subdivideECMesh->SetColumns((vtkDbiharStatic::SMC_CIRC / vtkDbiharStatic::EC_CIRC) * numSMCs);
	subdivideECMesh->Update();

	// vtkDbiharStatic::WritePolyData(subdivideECMesh->GetOutput(), (prefix + "quadMeshFullECc8000.vtp").c_str());
	vtkDbiharStatic::WritePolyData(subdivideECMesh->GetOutput(), "quadMeshFullECc8000.vtp");

	vtkSmartPointer<vtkSubdivideMesh> subdivideSMCMesh = vtkSmartPointer<vtkSubdivideMesh>::New();
	subdivideSMCMesh->SetInputData(fullMeshJoiner->GetOutput());
	subdivideSMCMesh->SetRows((vtkDbiharStatic::EC_AXIAL / vtkDbiharStatic::SMC_AXIAL) * numECs);
	subdivideSMCMesh->SetColumns(numSMCs);
	subdivideSMCMesh->Update();

	// vtkDbiharStatic::WritePolyData(subdivideSMCMesh->GetOutput(), (prefix + "quadMeshFullSMCc8000.vtp").c_str());
	vtkDbiharStatic::WritePolyData(subdivideSMCMesh->GetOutput(), "quadMeshFullSMCc8000.vtp");

#if 0
	int skipLength = 10;
	for (int i = 0; i < endPointIds->GetNumberOfIds(); i++)
	{
		vtkSmartPointer<vtkSkipSegmentFilter> skipSegmentFilter = vtkSmartPointer<vtkSkipSegmentFilter>::New();
		skipSegmentFilter->SetInputData(scalarRadiiToVectorsFilter->GetOutput());

		// Process only inlet on first iteration.
		skipSegmentFilter->SetInlet(i == 0);

		// Process only outlet on consecutive iterations.
		skipSegmentFilter->SetOutlet(i != 0);
		skipSegmentFilter->SetSkipSize(skipLength);
		skipSegmentFilter->SetPointId(endPointIds->GetId(i));
		skipSegmentFilter->SetNumberOfRadialQuads(numRadialQuads);
		skipSegmentFilter->Update();

		fullMeshJoiner->AddInputData(skipSegmentFilter->GetOutput());

		vtkSmartPointer<vtkEndCapFilter> endCapFilter = vtkSmartPointer<vtkEndCapFilter>::New();
		endCapFilter->SetInputData(skipSegmentFilter->GetOutput());
		endCapFilter->Update();

		vtkSmartPointer<vtkTriangleFilter> triangleCapFilter = vtkSmartPointer<vtkTriangleFilter>::New();
		triangleCapFilter->SetInputData(endCapFilter->GetOutput());
		triangleCapFilter->Update();

		// vtkDbiharStatic::WriteStlData(triangleCapFilter->GetOutput(), (prefix + SSTR("triCap_" << i << "_c8000.stl")).c_str());
		vtkDbiharStatic::WriteStlData(triangleCapFilter->GetOutput(), (SSTR("triCap_" << i << "_c8000.stl")).c_str());
	}

	fullMeshJoiner->Update();

	vtkSmartPointer<vtkTriangleFilter> triangleFullMeshFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFullMeshFilter->AddInputData(fullMeshJoiner->GetOutput());

	triangleFullMeshFilter->Update();
	// vtkDbiharStatic::WriteStlData(triangleFullMeshFilter->GetOutput(), (prefix + "triMeshFullc8000.stl").c_str());
	vtkDbiharStatic::WriteStlData(triangleFullMeshFilter->GetOutput(), "triMeshFullc8000.stl");
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}

