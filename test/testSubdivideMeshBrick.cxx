#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointSet.h>
#include <vtkSubdivideQuadFilter.h>
#include <vtkSubdivideQuadBrick.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkIdList.h>
#include <vtkAppendPolyData.h>

#include "wrapDbiharConfig.h"
#include "vtkDbiharStatic.h"
#include "vtkSubdivideMeshBrick.h"
#include "vtkJoinSmcBrickMesh.h"
#include <vtkTriangleFilter.h>

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLPolyDataReader> pointsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/quadMeshFullc4080.vtp").c_str());

	pointsReader->Update();

	vtkSmartPointer<vtkSubdivideMeshBrick> SubdivideMeshBrick = vtkSmartPointer<vtkSubdivideMeshBrick>::New();
	SubdivideMeshBrick->SetInputData(pointsReader->GetOutput());
	SubdivideMeshBrick->SetRows(52);
	SubdivideMeshBrick->SetColumns(4);
	SubdivideMeshBrick->SetAxialQuads(34);
	SubdivideMeshBrick->SetCircQuads(40);
	SubdivideMeshBrick->SetFlat(false);
	SubdivideMeshBrick->SetBranches(3);
	SubdivideMeshBrick->SetCellType(vtkDbiharStatic::SMC);
	SubdivideMeshBrick->Update();


	vtkDbiharStatic::WritePolyData(SubdivideMeshBrick->GetOutput(), "quadMeshFullc4080_SMCs.vtp");

	//vtkDbiharStatic::ShowPolyData(SubdivideMeshBrick->GetOutput());

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
