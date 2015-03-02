#include <map>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnsignedIntArray.h>
#include <vtkCellArray.h>

#include <vtkPointData.h>

#include <vtkXMLPolyDataWriter.h>
#include <vtkGenericDataObjectWriter.h>

#include <vtkXMLStructuredGridReader.h>

#include "showPolyData.h"
#include "vtkCentrelineToDbiharPatch.h"
#include "vtkCentrelinePartitioner.h"
#include "vtkScalarRadiiToVectorsFilter.h"

#include "vtkCentrelineData.h"
#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/227A_CentrelineResampled_4ECs.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(vesselCentreline);
	scalarRadiiToVectorsFilter->Update();

	vtkSmartPointer<vtkCentrelinePartitioner> centrelinePartitioner = vtkSmartPointer<vtkCentrelinePartitioner>::New();
	centrelinePartitioner->SetInputData(scalarRadiiToVectorsFilter->GetOutput());
	centrelinePartitioner->SetPartitionLength(50);
	vtkSmartPointer<vtkIdList> EndPoints = vtkSmartPointer<vtkIdList>::New();
	EndPoints->InsertNextId(20); // Remove the curve in inlet, causes folding in mesh otherwise.

	centrelinePartitioner->SetEndPoints(EndPoints);
	centrelinePartitioner->Update();

	vtkSmartPointer<vtkCentrelineToDbiharPatch> dbiharPatchFilter = vtkSmartPointer<vtkCentrelineToDbiharPatch>::New();
	dbiharPatchFilter->SetInputData(centrelinePartitioner->GetOutput());
	dbiharPatchFilter->SetNumberOfRadialQuads(28);
	dbiharPatchFilter->SetSpineId(0);
	dbiharPatchFilter->SetArchDerivScale(3.2);
	dbiharPatchFilter->SetEdgeDerivScale(4.0);
	dbiharPatchFilter->Update();

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputData(dbiharPatchFilter->GetOutput());
	writer->SetFileName("TEST_.vtk");
	writer->SetFileTypeToASCII();
	writer->Write();

	return EXIT_SUCCESS;
}

