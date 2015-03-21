#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkMath.h>

#include "vtkDbiharPatchFilter.h"
#include "vtkDbiharStatic.h"
#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	// Create the sides of a patch as 3D points.
	double x1 = -10.0;
	double x2 = 10.0;
	double y1 = -15.0;
	double y2 = 15.0;
	double z = 0;
	double arc = vtkMath::Pi();

	int cQuads = 26; // m = 25. Num quads should be even, to make sure m is odd.
	int yQuads = 30; // n = 29. Num quads should be even, to make sure n is odd.

	vtkIdType pIds = (cQuads + yQuads) * 2;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();

	// Insert boundary points. The boundary has four segments.
	// The coordinates of the points are calculated specific to the current boundary segment.
	double point[3] = {0.0};
	double radius = (x2 - x1) / 2;
	for(vtkIdType pId = 0; pId < pIds; pId++)
	{
		// Inserting points along the y = y1 boundary segment.
		if(pId < cQuads)
		{
			double dA = pId / (double)cQuads;
			if(dA < 0.5)
			{
				double angle = arc * dA;
				point[0] = -cos(angle) * radius;
				point[2] = sin(angle) * radius;
			}
			else
			{
				double angle = arc * dA - vtkMath::Pi() / 2.0;
				point[0] = sin(angle) * radius;
				point[2] = cos(angle) * radius;
			}
			point[1] = y1;
		}
		// Inserting points along the x = x2 boundary segment.
		else if(pId < cQuads + yQuads)
		{
			point[0] = x2;
			point[1] = y1 + (pId - cQuads) * ((y2 - y1) / yQuads);
			point[2] = z;
		}
		// Inserting points along the y = y2 boundary segment.
		else if(pId < cQuads * 2 + yQuads)
		{
			double dA = (pId - cQuads - yQuads) / (double)cQuads;
			if(dA < 0.5)
			{
				double angle = arc * dA;
				point[0] = cos(angle) * radius;
				point[2] = sin(angle) * radius;
			}
			else
			{
				double angle = arc * dA - vtkMath::Pi() / 2.0;
				point[0] = -sin(angle) * radius;
				point[2] = cos(angle) * radius;
			}
			point[1] = y2;
		}
		// Inserting points along the x = x1 boundary segment.
		else
		{
			point[0] = x1;
			point[1] = y2 - (pId - cQuads - yQuads - cQuads ) * ((y2 - y1) / yQuads);
			point[2] = z;
		}

		vtkIdType id = points->InsertNextPoint(point);
		// Sanity check.
		assert(id == pId);
		boundary->GetPointIds()->InsertNextId(pId);
	}
	boundary->GetPointIds()->InsertNextId(0);

	vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
	boundaries->InsertNextCell(boundary);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(points);
	inputPatch->SetLines(boundaries);

	vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
	derivatives->SetName(vtkDbiharPatchFilter::DERIV_ARR_NAME);
	derivatives->SetNumberOfComponents(3);

	double deriv[3] = {0.0};
	for(vtkIdType pId = 0; pId < pIds; pId++)
	{
		// Null them for the corner cases.
		deriv[0] = 0.0;
		deriv[1] = 0.0;
		deriv[2] = 0.0;

		// Inserting derivatives along the y = y1 boundary segment, skipping the corner case.
		if(pId < cQuads)
		{
			if(pId != 0)
			{
				deriv[0] = 0.0;
				deriv[1] = -10.0;
				deriv[2] = 0.0;
			}
		}
		// Inserting derivatives along the x = x2 boundary segment, skipping the corner case.
		else if(pId < cQuads + yQuads)
		{
			if(pId != cQuads)
			{
				deriv[0] = 0.0;
				deriv[1] = 0.0;
				deriv[2] = -60.0;
			}
		}
		// Inserting derivatives along the y = y2 boundary segment, skipping the corner case.
		else if(pId < cQuads * 2 + yQuads)
		{
			if(pId != cQuads + yQuads)
			{
				deriv[0] = 0.0;
				deriv[1] = 10.0;
				deriv[2] = 0.0;
			}
		}
		// Inserting derivatives along the x = x1 boundary segment, skipping the corner case.
		else
		{
			if(pId != cQuads * 2 + yQuads)
			{
				deriv[0] = 0.0;
				deriv[1] = 0.0;
				deriv[2] = -60.0;
			}
		}
		derivatives->InsertNextTuple(deriv);
	}
	derivatives->Print(std::cout);

	inputPatch->GetPointData()->SetVectors(derivatives);

	// showPolyData(inputPatch, NULL);

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	// Set the bounds of the UV space.
	patchFilter->SetA(0.0);
	patchFilter->SetB(2.0/3.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(vtkMath::Pi());
	// Set the boundary conditions.
	patchFilter->SetMQuads(cQuads);
	patchFilter->SetNQuads(yQuads);
	// Set solution method.
	patchFilter->SetIFlag(2);

	patchFilter->SetInputData(inputPatch);

	// patchFilter->Print(std::cout);

	patchFilter->Update();

	// patchFilter->Print(std::cout);

	vtkPolyData *outputPatch = patchFilter->GetOutput();

	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetDimensions(cQuads + 1, yQuads + 1, 1);
	structuredGrid->SetPoints(outputPatch->GetPoints());

	vtkDbiharStatic::ShowPolyDataWithGrid(inputPatch, structuredGrid);

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
