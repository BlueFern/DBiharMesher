#include <stdlib.h>

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>

#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>

#include "vtkDbiharPatchFilter.h"

void showPolyData(vtkPolyData *polyData);

int main(int argc, char* argv[]) {

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();
	patchFilter->SetITCG(4);

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	// Create the sides of a patch as 3D points.
	double x = -10.0;
	double width = 20.0;
	double y = -15.0;
	double height = 30;

	int xQuads = 20; // m = 19. Num quads should be even, to make sure m is odd.
	int yQuads = 30; // n = 29. Num quads should be even, to make sure n is odd.

	//points->SetNumberOfPoints((xQuads + 1) * 2 + (yQuads + 1) * 2);

	vtkSmartPointer<vtkPolyLine> lowBorder = vtkSmartPointer<vtkPolyLine>::New();
	vtkSmartPointer<vtkPolyLine> highBorder = vtkSmartPointer<vtkPolyLine>::New();

	// Inserting points for lines y and y + height.
	double dX = width / xQuads;
	for(int i = 0; i <= xQuads; i++)
	{
		double pX = x + i * dX;
		vtkIdType pId = points->InsertNextPoint(pX, y, 0.0);
		lowBorder->GetPointIds()->InsertNextId(pId);
	}

	for(int i = 0; i <= xQuads; i++)
	{
		double pX = x + i * dX;
		vtkIdType pId = points->InsertNextPoint(pX, y + height, 0.0);
		highBorder->GetPointIds()->InsertNextId(pId);
	}

	// Should be 42 points.
	std::cout << points->GetNumberOfPoints() << std::endl;

	vtkSmartPointer<vtkPolyLine> leftBorder = vtkSmartPointer<vtkPolyLine>::New();
	vtkSmartPointer<vtkPolyLine> rightBorder = vtkSmartPointer<vtkPolyLine>::New();

	// Inserting points for lines x and x + width.
	double dY = height / yQuads;
	for(int i = 0; i <= yQuads; i++)
	{
		double pY = y + i * dY;
		vtkIdType pId = points->InsertNextPoint(x, pY, 0.0);
		leftBorder->GetPointIds()->InsertNextId(pId);
	}

	for(int i = 0; i <= yQuads; i++)
	{
		double pY = y + i * dY;
		vtkIdType pId = points->InsertNextPoint(x + width, pY, 0.0);
		rightBorder->GetPointIds()->InsertNextId(pId);
	}

	// Should be 42 + 62 points.
	std::cout << points->GetNumberOfPoints() << std::endl;

	std::cout << lowBorder->GetNumberOfPoints() << ", " << highBorder->GetNumberOfPoints() << ", " << leftBorder->GetNumberOfPoints() << ", " << rightBorder->GetNumberOfPoints() << std::endl;

	vtkSmartPointer<vtkCellArray> borders = vtkSmartPointer<vtkCellArray>::New();

	borders->InsertNextCell(lowBorder);
	borders->InsertNextCell(highBorder);
	borders->InsertNextCell(leftBorder);
	borders->InsertNextCell(rightBorder);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(points);
	inputPatch->SetLines(borders);

	vtkIdList *pts = vtkIdList::New();
	inputPatch->GetLines()->GetCell(0, pts);
	pts->Print(std::cout);

	showPolyData(inputPatch);

	patchFilter->SetInputData(inputPatch);

	patchFilter->Print(std::cout);

	patchFilter->Update();

	std::cout << "Output points: " << patchFilter->GetOutput()->GetNumberOfPoints() << std::endl;

	showPolyData(patchFilter->GetOutput());

	patchFilter->Print(std::cout);

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}

void showPolyData(vtkPolyData *polyData)
{
	vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
	for(int pId = 0; pId < polyData->GetNumberOfPoints(); pId++)
	{
		idList->InsertNextId(pId);
	}

	vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
	verts->InsertNextCell(idList);

	polyData->SetVerts(verts);
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(polyData);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetPointSize(1);

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkInteractorStyleSwitch *iss = vtkInteractorStyleSwitch::SafeDownCast(renderWindowInteractor->GetInteractorStyle());
	iss->SetCurrentStyleToTrackballCamera();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderer->AddActor(actor);

	renderWindow->Render();
	renderWindowInteractor->Start();

	renderWindow->Finalize();
}
