#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointSet.h>
#include <vtkSubdivideQuadFilter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkIdList.h>

#include "vtkSubdivideMeshDynamic.h"
#include "wrapDbiharConfig.h"
#include "vtkDbiharStatic.h"


int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLPolyDataReader> pointsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/testSubdivideMeshDynamic0.vtp").c_str());

	pointsReader->Update();

	vtkSmartPointer<vtkSubdivideMeshDynamic> subdivideMeshDynamic = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	subdivideMeshDynamic->SetInputData(pointsReader->GetOutput());
	subdivideMeshDynamic->SetHeight(3);
	subdivideMeshDynamic->SetLength(3);
	subdivideMeshDynamic->Print(std::cout);
	subdivideMeshDynamic->Update();

	vtkDbiharStatic::ShowPolyData(subdivideMeshDynamic->GetOutput());

#if 0
	vtkDbiharStatic::WritePolyData(subdivideMeshDynamic->GetOutput(), "subdivideMeshDynamicTest.vtp");
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
