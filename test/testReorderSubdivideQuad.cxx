#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkIdList.h>

#include "vtkReorderSubdivideQuad.h"
#include "wrapDbiharConfig.h"
#include "vtkDbiharStatic.h"


int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLPolyDataReader> pointsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/testSubdivideQuadBrick.vtp").c_str());

	pointsReader->Update();
	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	pd->DeepCopy(pointsReader->GetOutput());
	vtkSmartPointer<vtkReorderSubdivideQuad> reorderSubdivideQuad = vtkSmartPointer<vtkReorderSubdivideQuad>::New();
	vtkSmartPointer<vtkPolyData> pointsData = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkIdList> pointsList = vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	pd->GetCellPoints(0,pointsList);
	points->SetNumberOfPoints(4);
	points->SetPoint(0,pd->GetPoint(pointsList->GetId(0)));
	points->SetPoint(1,pd->GetPoint(pointsList->GetId(1)));
	points->SetPoint(2,pd->GetPoint(pointsList->GetId(2)));
	points->SetPoint(3,pd->GetPoint(pointsList->GetId(3)));
	pointsData->SetPoints(points);
	reorderSubdivideQuad->SetInputData(pointsData);
	reorderSubdivideQuad->SetRows(20);
	reorderSubdivideQuad->SetColumns(4);

	reorderSubdivideQuad->SetRotations(1); // CLOCKWISE rotations - TODO: write some doc pls
	reorderSubdivideQuad->Update();

	vtkDbiharStatic::ShowPolyData(reorderSubdivideQuad->GetOutput());
#if 1
	vtkDbiharStatic::WritePolyData(reorderSubdivideQuad->GetOutput(), "reorder_out.vtp");
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
