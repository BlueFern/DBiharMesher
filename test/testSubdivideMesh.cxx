#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointSet.h>
#include <vtkSubdivideQuadFilter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkIdList.h>
#include <vtkAppendPolyData.h>

#include "wrapDbiharConfig.h"
#include "vtkDbiharStatic.h"
#include "vtkSubdivideMesh.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLPolyDataReader> pointsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/output.vtp").c_str());

	pointsReader->Update();



	vtkSmartPointer<vtkSubdivideMesh> subdivideMesh = vtkSmartPointer<vtkSubdivideMesh>::New();
	subdivideMesh->SetInputData(pointsReader->GetOutput());
	subdivideMesh->SetRows(52);
	subdivideMesh->SetColumns(4);
	subdivideMesh->Update();
	vtkDbiharStatic::WritePolyData(subdivideMesh->GetOutput(), "quadMeshFullSMCc4080.vtp");

	vtkDbiharStatic::ShowPolyData(subdivideMesh->GetOutput());

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
