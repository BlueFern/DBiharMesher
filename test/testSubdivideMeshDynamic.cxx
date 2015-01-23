#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointSet.h>
#include <vtkSubdivideQuadFilter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkIdList.h>

#include "vtkSubdivideMeshDynamic.h"
#include "wrapDbiharConfig.h"
#include "showPolyData.h"


int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLPolyDataReader> pointsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/mesh0.vtp").c_str());

	pointsReader->Update();

	vtkSmartPointer<vtkSubdivideMeshDynamic> subdivideMeshDynamic = vtkSmartPointer<vtkSubdivideMeshDynamic>::New();
	subdivideMeshDynamic->SetInputData(pointsReader->GetOutput());
	subdivideMeshDynamic->SetHeight(0.05);
	subdivideMeshDynamic->SetLength(1);
	subdivideMeshDynamic->Print(std::cout);
	subdivideMeshDynamic->Update();

#if 1
	vtkSmartPointer<vtkXMLPolyDataWriter> meshWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	meshWriter->SetInputData(subdivideMeshDynamic->GetOutput());
	meshWriter->SetFileName("dynamicTinyMesh0.vtp");
	meshWriter->Update();


#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
