#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vktPolyData.h>

#include "vtkDbiharPatchFilter.h"

#include "showPolyData.h"

int main(int argc, char* argv[]) {

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	patchFilter->SetA(0.0);
	patchFilter->SetB(2.0/3.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(1.0);

	patchFilter->SetIFlag(2);

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	// Create the sides of a patch as 3D points.
	double x = -10.0;
	double width = 20.0;
	double y = -15.0;
	double height = 30;

	int xQuads = 20; // m = 19. Num quads should be even, to make sure m is odd.
	int yQuads = 30; // n = 29. Num quads should be even, to make sure n is odd.

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

	vtkSmartPointer<vtkCellArray> borders = vtkSmartPointer<vtkCellArray>::New();

	borders->InsertNextCell(lowBorder);
	borders->InsertNextCell(highBorder);
	borders->InsertNextCell(leftBorder);
	borders->InsertNextCell(rightBorder);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(points);
	inputPatch->SetLines(borders);

	patchFilter->SetInputData(inputPatch);

	patchFilter->Print(std::cout);

	patchFilter->Update();

	patchFilter->Print(std::cout);

	vtkPolyData *outputPatch = patchFilter->GetOutput();

	showPolyData(inputPatch, outputPatch);

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
