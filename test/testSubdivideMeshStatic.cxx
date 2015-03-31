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
#include "vtkSubdivideMeshStatic.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLPolyDataReader> pointsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/mesh0.vtp").c_str());

	pointsReader->Update();

	vtkDbiharStatic::ShowPolyData(pointsReader->GetOutput());
	vtkSmartPointer<vtkSubdivideMeshStatic> subdivideMeshStatic = vtkSmartPointer<vtkSubdivideMeshStatic>::New();
	subdivideMeshStatic->SetInputData(pointsReader->GetOutput());
	subdivideMeshStatic->SetColumns(2);
	subdivideMeshStatic->SetRows(1);
	subdivideMeshStatic->Print(std::cout);
	subdivideMeshStatic->Update();

	vtkDbiharStatic::ShowPolyData(subdivideMeshStatic->GetOutput());

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
