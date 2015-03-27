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

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLPolyDataReader> pointsReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	pointsReader->SetFileName((std::string(TEST_DATA_DIR) + "/mesh0.vtp").c_str());

	pointsReader->Update();
	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	pd->DeepCopy(pointsReader->GetOutput());

	vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
	vtkSmartPointer<vtkIdList> pointsList = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	pointsList->SetNumberOfIds(4);
	points->SetNumberOfPoints(4);

	for (int i = 0; i < pd->GetNumberOfCells(); i++)
	{
		vtkSmartPointer<vtkSubdivideQuadFilter> subdivideQuadFilter = vtkSmartPointer<vtkSubdivideQuadFilter>::New();
		vtkSmartPointer<vtkPolyData> pointsData = vtkSmartPointer<vtkPolyData>::New();
		pd->GetCellPoints(i,pointsList);
		points->SetPoint(0,pd->GetPoint(pointsList->GetId(0)));
		points->SetPoint(1,pd->GetPoint(pointsList->GetId(1)));
		points->SetPoint(2,pd->GetPoint(pointsList->GetId(2)));
		points->SetPoint(3,pd->GetPoint(pointsList->GetId(3)));
		pointsData->SetPoints(points);
		subdivideQuadFilter->SetInputData(pointsData);
		subdivideQuadFilter->SetColumns(4);
		subdivideQuadFilter->SetRows(52);
		subdivideQuadFilter->Update();
		appendPolyData->AddInputData(subdivideQuadFilter->GetOutput());
	}
	appendPolyData->Update();

#if 1
	vtkDbiharStatic::WritePolyData(appendPolyData->GetOutput(), "subdivideMeshStaticTest.vtp");
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
