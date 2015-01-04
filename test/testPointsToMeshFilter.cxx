#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "vtkPointsToMeshFilter.h"
#include "showPolyData.h"


int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLPolyDataReader> pointsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	pointsReader->SetFileName("testPoints.vtp");

	vtkSmartPointer<vtkUnsignedIntArray> dimensions = vtkSmartPointer<vtkUnsignedIntArray>::New();
	dimensions->InsertNextValue(28);
	dimensions->InsertNextValue(44);
	dimensions->InsertNextValue(44);
	dimensions->InsertNextValue(44);
	dimensions->InsertNextValue(44);
	dimensions->InsertNextValue(44);
	dimensions->InsertNextValue(44);

	vtkSmartPointer<vtkPointsToMeshFilter> pointsToMeshFilter = vtkSmartPointer<vtkPointsToMeshFilter>::New();
	//pointsToMeshFilter->DebugOn();
	pointsToMeshFilter->SetDimensions(dimensions);
	pointsToMeshFilter->SetInputConnection(pointsReader->GetOutputPort());

	pointsToMeshFilter->Update();

	showPolyData1(pointsToMeshFilter->GetOutput(), 1.0);

	std::cout << "Number of points in the output: " << pointsToMeshFilter->GetOutput()->GetNumberOfPoints() << std::endl;

#if 0
	vtkSmartPointer<vtkXMLPolyDataWriter> meshWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	meshWriter->SetInputConnection(meshWriter->GetOutputPort());
	meshWriter->SetFileName("outputMesh.vtp");
	meshWriter->Update();
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
