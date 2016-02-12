#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkGenericDataObjectReader.h>

#include "wrapDbiharConfig.h"
#include "vtkDbiharPatchSmooth.h"
#include "vtkDbiharStatic.h"
#include "vtkRescaleUnits.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;


	vtkSmartPointer<vtkXMLPolyDataReader> pointsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/quadMeshFullc4080.vtp").c_str());

	pointsReader->Update();
	vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::New();
	mesh->DeepCopy(pointsReader->GetOutput());

	vtkSmartPointer<vtkDbiharPatchSmooth> dbiharPatchSmooth = vtkSmartPointer<vtkDbiharPatchSmooth>::New();
	dbiharPatchSmooth->SetInputData(mesh);
	dbiharPatchSmooth->SetNumRadialQuads(20);

	dbiharPatchSmooth->Update();
	vtkDbiharStatic::WritePolyData(dbiharPatchSmooth->GetOutput(), "new_Bifurcation_Model.vtp");
	vtkDbiharStatic::ShowPolyData(dbiharPatchSmooth->GetOutput());


	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
