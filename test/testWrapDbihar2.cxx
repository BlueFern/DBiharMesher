#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkMath.h>

#include "vtkDbiharPatchFilter.h"

#include "showPolyData.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	patchFilter->SetA(0.0);
	patchFilter->SetB(2.0/3.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(vtkMath::Pi());

	patchFilter->SetIFlag(2);

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	// Create the sides of a patch as 3D points.
	double x = -10.0;
	double arc = vtkMath::Pi();
	double rad = 10.0;
	double y = -15.0;
	double height = 30;

	int cQuads = 26; // m = 25. Num quads should be even, to make sure m is odd.
	int yQuads = 30; // n = 29. Num quads should be even, to make sure n is odd.

	vtkSmartPointer<vtkPolyLine> lowBorder = vtkSmartPointer<vtkPolyLine>::New();
	vtkSmartPointer<vtkPolyLine> highBorder = vtkSmartPointer<vtkPolyLine>::New();

	// Inserting points for lines y and y + height.
	for(int i = 0; i <= cQuads; i++)
	{
		double dA = i / (double)cQuads;
		if(dA < 0.5)
		{
			double angle = arc * dA;
			double px = -(cos(angle) * rad);
			double pz = sin(angle) * rad;
			lowBorder->GetPointIds()->InsertNextId(points->InsertNextPoint(px, y, pz));
			highBorder->GetPointIds()->InsertNextId(points->InsertNextPoint(px, y + height, pz));
		}
		else
		{
			double angle = arc * dA - vtkMath::Pi() / 2.0;
			double px = sin(angle) * rad;
			double pz = cos(angle) * rad;
			lowBorder->GetPointIds()->InsertNextId(points->InsertNextPoint(px, y, pz));
			highBorder->GetPointIds()->InsertNextId(points->InsertNextPoint(px, y + height, pz));
		}
	}

	vtkSmartPointer<vtkPolyLine> leftBorder = vtkSmartPointer<vtkPolyLine>::New();
	vtkSmartPointer<vtkPolyLine> rightBorder = vtkSmartPointer<vtkPolyLine>::New();

	// Inserting points for lines x and x + width.
	double dY = height / yQuads;
	for(int i = 0; i <= yQuads; i++)
	{
		double pY = y + i * dY;
		leftBorder->GetPointIds()->InsertNextId(points->InsertNextPoint(x, pY, 0.0));
		rightBorder->GetPointIds()->InsertNextId(points->InsertNextPoint(x + rad * 2, pY, 0.0));
	}

	vtkSmartPointer<vtkCellArray> borders = vtkSmartPointer<vtkCellArray>::New();

	borders->InsertNextCell(lowBorder);
	borders->InsertNextCell(highBorder);
	borders->InsertNextCell(leftBorder);
	borders->InsertNextCell(rightBorder);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(points);
	inputPatch->SetLines(borders);

	// showPolyData(inputPatch, NULL);

	patchFilter->SetInputData(inputPatch);

	patchFilter->Print(std::cout);

	patchFilter->Update();

	patchFilter->Print(std::cout);

	vtkPolyData *outputPatch = patchFilter->GetOutput();

	showPolyData(inputPatch, outputPatch);

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
