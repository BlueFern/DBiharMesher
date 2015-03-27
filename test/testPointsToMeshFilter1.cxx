#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>

#include "wrapDbiharConfig.h"
#include "vtkPointsToMeshFilter.h"
#include "vtkDbiharStatic.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLPolyDataReader> pointsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/dbhPatchPoints0.vtp").c_str());
	//pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/dbhPatchPoints1.vtp").c_str());
	//pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/dbhPatchPoints2.vtp").c_str());

	vtkSmartPointer<vtkUnsignedIntArray> dimensions = vtkSmartPointer<vtkUnsignedIntArray>::New();

	// Straight segment.
	dimensions->InsertNextValue(32); // Use dbhPatchPoints0.vtp
	dimensions->InsertNextValue(70);

	// Even branch length bifurcations.
	//dimensions->InsertNextValue(28);
	//dimensions->InsertNextValue(44);
	//dimensions->InsertNextValue(44); // Use dbhPatchPoints1.vtp
	//dimensions->InsertNextValue(44);

	// Different length branches bifurcation.
	//dimensions->InsertNextValue(28);
	//dimensions->InsertNextValue(70);
	//dimensions->InsertNextValue(24); // Use dbhPatchPoints2.vtp
	//dimensions->InsertNextValue(44);

	vtkSmartPointer<vtkPointsToMeshFilter> pointsToMeshFilter = vtkSmartPointer<vtkPointsToMeshFilter>::New();
	//pointsToMeshFilter->DebugOn();
	pointsToMeshFilter->SetDimensions(dimensions);
	pointsToMeshFilter->SetInputConnection(pointsReader->GetOutputPort());
	pointsToMeshFilter->SetShowProgress(true);
	pointsToMeshFilter->Update();
	pointsToMeshFilter->Print(std::cout);

	//showPolyData1(pointsToMeshFilter->GetOutput(), 1.0);

	std::cout << "Number of points in the output: " << pointsToMeshFilter->GetOutput()->GetNumberOfPoints() << std::endl;

	vtkDbiharStatic::WritePolyData(pointsToMeshFilter->GetOutput(), "meshTestOutput.vtp");

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
