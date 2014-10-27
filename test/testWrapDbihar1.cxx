#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkStructuredGrid.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>

#include "vtkDbiharPatchFilter.h"
#include "showPolyData.h"

#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	// Create the sides of a patch as 3D points.
	double x1 = -10.0;
	double x2 = 10.0;
	double y1 = -15.0;
	double y2 = 15.0;
	double z = 0.0;

	int xQuads = 20; // m = 19. Num quads should be even, to make sure m is odd.
	int yQuads = 30; // n = 29. Num quads should be even, to make sure n is odd.

	vtkIdType pIds = (xQuads + yQuads) * 2;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();

	// Insert boundary points. The boundary has four segments.
	// The coordinates of the points are calculated specific to the current boundary segment.
	double point[3] = {0.0, 0.0, 0.0};
	for(vtkIdType pId = 0; pId < pIds; pId++)
	{
		// Inserting along the y = y1 boundary segment.
		if(pId < xQuads)
		{
			point[0] = x1 + pId * ((x2 - x1) / xQuads);
			point[1] = y1;
		}
		// Inserting along the x = x2 boundary segment.
		else if(pId < xQuads + yQuads)
		{
			point[0] = x2;
			point[1] = y1 + (pId - xQuads) * ((y2 - y1) / yQuads);
		}
		// Inserting along the y = y2 boundary segment.
		else if(pId < xQuads * 2 + yQuads)
		{
			point[0] = x2 - (pId - xQuads - yQuads) * ((x2 - x1) / xQuads);
			point[1] = y2;
		}
		// Inserting along the x = x1 boundary segment.
		else
		{
			point[0] = x1;
			point[1] = y2 - (pId - xQuads - yQuads - xQuads ) * ((y2 - y1) / yQuads);
		}

		// Z coordinate does not change in this case.
		point[2] = z;

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

	// showPolyData(inputPatch, NULL);

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	// Set the bounds of the UV space.
	patchFilter->SetA(0.0);
	patchFilter->SetB(2.0/3.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(1.0);
	// Set the boundary conditions.
	patchFilter->SetMQuads(xQuads);
	patchFilter->SetNQuads(yQuads);
	// Set solution method.
	patchFilter->SetIFlag(2);

	patchFilter->SetInputData(inputPatch);

	// patchFilter->Print(std::cout);

	patchFilter->Update();

	//  patchFilter->Print(std::cout);

	vtkPolyData *outputPatch = patchFilter->GetOutput();

	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetDimensions(xQuads + 1, yQuads + 1, 1);
	structuredGrid->SetPoints(outputPatch->GetPoints());

	showPolyData(inputPatch, structuredGrid);

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
