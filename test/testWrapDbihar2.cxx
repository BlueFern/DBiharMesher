#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkMath.h>

#include "vtkDbiharPatchFilter.h"
#include "vtkDbiharStatic.h"
#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	// Create the sides of a patch as 3D points.
	double radius = 10;
	double length = 30;

	double a = 1.0; // vtkMath::Pi();

	int cQuads = 26; // m = 25. Num quads should be even, to make sure m is odd.
	int xQuads = 30; // n = 29. Num quads should be even, to make sure n is odd.

	vtkIdType pIds = (cQuads + xQuads) * 2;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();

	// Insert boundary points. The boundary has four segments.
	// The coordinates of the points are calculated specific to the current boundary segment.
	double point[3] = {0.0};

	vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
	derivatives->SetName(vtkDbiharPatchFilter::DERIV_ARR_NAME);
	derivatives->SetNumberOfComponents(3);

	double S1 = 12; // Horizontal boundaries. (Along the arch x = +/- 15)
	double S2 = 39.5; // Vertical boundaries. (Along straight Line y = +/- 10)

	double deriv[3] = {0.0};

	for(vtkIdType pId = 0; pId < pIds; pId++)
	{
		// Null them for the corner cases.
		deriv[0] = 0.0;
		deriv[1] = 0.0;
		deriv[2] = 0.0;

		// Inserting points and derivatives along the u = 0 boundary segment.
		if(pId < cQuads)
		{
			double V = vtkMath::Pi() * (double)pId / (double)cQuads;
			point[0] = 0;
			point[1] = radius * cos(a * V);
			point[2] = radius * sin(a * V);

			if(pId != 0)
			{
				deriv[0] = -S1;
				deriv[1] = 0.0;
				deriv[2] = 0.0;
			}
		}
		// Inserting points along the v = 1 boundary segment.
		else if(pId < cQuads + xQuads)
		{
			point[0] = (pId - cQuads) * (length / xQuads);
			point[1] = -radius;
			point[2] = 0;

			if(pId != cQuads)
			{
				deriv[0] = 0.0;
				deriv[1] = 0.0;
				deriv[2] = -S2;
			}
		}
		// Inserting points along the u = 1 boundary segment.
		else if(pId < cQuads * 2 + xQuads)
		{
			double V = vtkMath::Pi() * (double)(pId - cQuads - xQuads) / (double)cQuads;

			point[0] = length;
			point[1] = -radius * cos(a * V);
			point[2] = radius * sin(a * V);

			if(pId != cQuads + xQuads)
			{
				deriv[0] = S1;
				deriv[1] = 0.0;
				deriv[2] = 0.0;
			}
		}
		// Inserting points along the v = 0 boundary segment.
		else
		{
			point[0] = length - (pId - (2*cQuads + xQuads)) * (length / xQuads);
			point[1] = radius;
			point[2] = 0;

			if(pId != cQuads * 2 + xQuads)
			{
				deriv[0] = 0.0;
				deriv[1] = 0.0;
				deriv[2] = -S2;
			}
		}

		vtkIdType id = points->InsertNextPoint(point);

		// Sanity check.
		assert(id == pId);
		boundary->GetPointIds()->InsertNextId(pId);

		derivatives->InsertNextTuple(deriv);
	}
	boundary->GetPointIds()->InsertNextId(0);

	vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
	boundaries->InsertNextCell(boundary);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(points);
	inputPatch->SetLines(boundaries);

	inputPatch->GetPointData()->SetVectors(derivatives);

	vtkDbiharStatic::ShowPolyData(inputPatch);

	vtkDbiharStatic::WritePolyData(inputPatch, std::string(argv[0]) + "_inputPatch.vtp");

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	// Set the bounds of the UV space.
	patchFilter->SetA(0.0);
	patchFilter->SetB(1);
	patchFilter->SetC(0.0);
	patchFilter->SetD(vtkMath::Pi());
	// Set the number of quads.
	patchFilter->SetMQuads(cQuads);
	patchFilter->SetNQuads(xQuads);
	// Set the boundary conditions.
	patchFilter->SetInputData(inputPatch);
	// Set solution method.
	patchFilter->SetIFlag(2);

	patchFilter->Update();

	vtkPolyData *outputPatch = patchFilter->GetOutput();

	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetDimensions(cQuads + 1, xQuads + 1, 1);
	structuredGrid->SetPoints(outputPatch->GetPoints());

	vtkDbiharStatic::ShowPolyDataWithGrid(inputPatch, structuredGrid);

	// Convert vtkStructuredGrid object to vtkPolydata.
	vtkSmartPointer<vtkStructuredGridGeometryFilter> gridGeometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	gridGeometryFilter->SetInputData(structuredGrid);
	gridGeometryFilter->Update();

	vtkDbiharStatic::WritePolyData(gridGeometryFilter->GetOutput(), std::string(argv[0]) + "_outputSurface.vtp");

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
