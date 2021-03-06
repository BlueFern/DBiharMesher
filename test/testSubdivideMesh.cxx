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
	pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/testSubdivideMesh.vtp").c_str());

	pointsReader->Update();

	vtkDbiharStatic::WritePolyData(pointsReader->GetOutput(), "t1.vtp");

	vtkSmartPointer<vtkSubdivideMesh> subdivideMesh = vtkSmartPointer<vtkSubdivideMesh>::New();
	subdivideMesh->SetInputData(pointsReader->GetOutput());
	subdivideMesh->SetColumns(5);
	subdivideMesh->SetRows(1);
	subdivideMesh->Print(std::cout);
	subdivideMesh->Update();

	vtkDbiharStatic::ShowPolyData(subdivideMesh->GetOutput());

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
