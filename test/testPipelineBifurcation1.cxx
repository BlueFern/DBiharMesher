#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSTLWriter.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkGenericCell.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkAppendPoints.h>
#include <vtkAppendPolyData.h>
#include <vtkTriangleFilter.h>

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
#include "showPolyData.h"

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

	vtkSmartPointer<vtkCentrelinePartitioner> centrelinePartitioner = vtkSmartPointer<vtkCentrelinePartitioner>::New();
	centrelinePartitioner->SetInputData(scalarRadiiToVectorsFilter->GetOutput());
	centrelinePartitioner->SetPartitionLength(30);
	centrelinePartitioner->Update();

#if 0
	vtkSmartPointer<vtkXMLPolyDataWriter> tmpWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	tmpWriter->SetFileName("tmpCentreline.vtp");
	tmpWriter->SetInputData(centrelinePartitioner->GetOutput());
	tmpWriter->Update();

	std::cout << "Exiting " << __FUNCTION__ << " in " << __FILE__ << ":" << __LINE__ << std::endl;
	exit(EXIT_FAILURE);
#endif

	vtkPolyData *partitionedCentreline = centrelinePartitioner->GetOutput();

	//double lengths[partitionedCentreline->GetLines()->GetNumberOfCells()] = {0.0};
	std::vector<double> lengths(partitionedCentreline->GetLines()->GetNumberOfCells());

	vtkSmartPointer<vtkAppendPoints> appendPoints = vtkSmartPointer<vtkAppendPoints>::New();

	// Working with centreline partitions 3, 4, 5.
	for (int i = 9; i < 12; i++)
	{
		vtkSmartPointer<vtkCentrelineToDbiharPatch> dbiharPatchFilter = vtkSmartPointer<vtkCentrelineToDbiharPatch>::New();
		dbiharPatchFilter->SetInputData(partitionedCentreline);
		dbiharPatchFilter->SetNumberOfRadialQuads(28);
		dbiharPatchFilter->SetSpineId(i);
		dbiharPatchFilter->SetArchDerivScale(3.2);
		dbiharPatchFilter->SetEdgeDerivScale(4.0);
		dbiharPatchFilter->Update();

		// WARNING: This statement is overwriting memory, because the size of lengths is only 3. Here we are writing past 3.
		lengths[i] = partitionedCentreline->GetCell(i)->GetNumberOfPoints();
		appendPoints->AddInputData(dbiharPatchFilter->GetOutput());
	}

	appendPoints->Update();

	vtkSmartPointer<vtkUnsignedIntArray> dimensions = vtkSmartPointer<vtkUnsignedIntArray>::New();
	dimensions->InsertNextValue(28);

	// Solving simultaneous equations for each branch sections length.
	dimensions->InsertNextValue((lengths[9] - lengths[10] + lengths[11]) / 2);
	dimensions->InsertNextValue((lengths[9] + lengths[10] - lengths[11]) / 2);
	dimensions->InsertNextValue((-lengths[9] + lengths[10] + lengths[11]) / 2);

	vtkSmartPointer<vtkPointsToMeshFilter> pointsToMeshFilter = vtkSmartPointer<vtkPointsToMeshFilter>::New();
	pointsToMeshFilter->SetInputData(appendPoints->GetOutput());
	pointsToMeshFilter->SetDimensions(dimensions);
	pointsToMeshFilter->Update();

	showPolyData1(pointsToMeshFilter->GetOutput());

	vtkSmartPointer<vtkXMLPolyDataWriter> quadMeshWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	quadMeshWriter->SetInputData(pointsToMeshFilter->GetOutput());
	quadMeshWriter->SetFileName("quadMeshBifurcation1.vtp");
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
	ECMeshWriter->SetFileName("ECquadMeshBifurcation1.vtp");
	ECMeshWriter->Write();

	vtkSmartPointer<vtkSubdivideMeshDynamic> dynamicSMCMesher = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	dynamicSMCMesher->SetInputData(pointsToMeshFilter->GetOutput());
	dynamicSMCMesher->SetHeight(vtkDbiharStatic::SMC_CIRC);
	dynamicSMCMesher->SetLength(vtkDbiharStatic::SMC_AXIAL);
	dynamicSMCMesher->Update();

	vtkSmartPointer<vtkXMLPolyDataWriter> SMCMeshWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	SMCMeshWriter->SetInputData(dynamicSMCMesher->GetOutput());
	SMCMeshWriter->SetFileName("SMCquadMeshBifurcation1.vtp");
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
	vtkSmartPointer<vtkIdList> endPointIds = vtkSmartPointer<vtkIdList>::New();
	for (int i = 3; i < 6; i++)
	{
		endPointIds->InsertNextId(partitionedCentreline->GetCell(i)->GetPointId(0));
	}

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
	triMeshWriter->SetFileName("triMeshWithCapsBifurcation1.stl");
	triMeshWriter->Write();

	return EXIT_SUCCESS;
}
