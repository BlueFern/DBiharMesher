/*
 * Tube/Cylinder Patch Generator
 */

#include <stdlib.h>

#include <sstream>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkTriangleFilter.h>

#include "vtkDbiharPatchFilter.h"
#include "vtkDbiharStatic.h"

#include "wrapDbiharConfig.h"

#define SSTR( x ) dynamic_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()

int main(int argc, char* argv[])
{

	std::cout << "Starting " << __FILE__ << std::endl;

	int xQuads = 24; // m = 23. Num quads should be even, to make sure m is odd.
	int yQuads = 30; // n = 29. Num quads should be even, to make sure n is odd.
	int numPoints = (xQuads + yQuads) * 2;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(numPoints);

	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();
	boundary->GetPointIds()->SetNumberOfIds(numPoints);

	// INSERT POINTS HERE.
	// Insert boundary points. The boundary has four segments.
	double point[3] = {0.0};
	double angleInc = ((vtkMath::Pi()) * 2) / (xQuads);
	double radius = 5;
	double length = 25;

	for(vtkIdType pId = 0; pId < numPoints; pId++ )
	{

		// Segment 1: first Circle (y = 0)
		if (pId < xQuads)
		{
			point[0] = cos(pId * angleInc) * radius;
			point[1] = - (length / 2);
			point[2] = sin(pId * angleInc) * radius;
		}

		// Segment 2: Straight Line up at (x = Radius) and (z = 0)
		else if (pId < (xQuads + yQuads))
		{
			point[0] = radius;
			point[1] = (- length /2) + (pId - xQuads) * (length / yQuads);
			point[2] = 0.0;
		}

		// Segment 3: second Circle (y = length)
		else if (pId < (2 * xQuads + yQuads))
		{
			point[0] = cos((pId - (xQuads + yQuads)) * angleInc) * radius;
			point[1] = length / 2;
			point[2] = -sin((pId - (xQuads + yQuads)) * angleInc) * radius;
		}

		// Segment 4: Straight Line down at (x = Radius) and (z = 0)
		else
		{
			point[0] = radius;
			point[1] = (length / 2) - ((pId - (2 * xQuads + yQuads)) * (length / yQuads));
			point[2] = 0.0;
		}


		points->InsertPoint(pId, point);
		boundary->GetPointIds()->InsertId(pId, pId);
	}

	boundary->GetPointIds()->InsertId(numPoints, 0);

	for(vtkIdType pId = 0; pId < boundary->GetPointIds()->GetNumberOfIds(); pId++)
	{
		std::cout << boundary->GetPointIds()->GetId(pId) << " ";
	}

	vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
	boundaries->InsertNextCell(boundary);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(points);
	vtkDbiharStatic::WritePolyData(inputPatch, std::string(argv[0]) + "_points.vtp");
	inputPatch->SetLines(boundaries);

	// INSERT POINTS DONE.

	// INSERT DERIVATIVES HERE.

	vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
	derivatives->SetName(vtkDbiharPatchFilter::DERIV_ARR_NAME);
	derivatives->SetNumberOfComponents(3);

	double deriv[3] = {0.0};
	double dY = 63;
	double dZ = 16;
	double B = 1.9;
	double D = 0.4;

	for (vtkIdType pId_d = 0; pId_d < numPoints; pId_d++)
	{

		// Derivatives on Edge-Points = {0.0, 0.0, 0.0}
		if ((pId_d == 0) || (pId_d == xQuads) || (pId_d == (xQuads + yQuads)) || (pId_d == (2 * xQuads + yQuads)))
		{
			deriv[0] = 0.0;
			deriv[1] = 0.0;
			deriv[2] = 0.0;
		}

		// Derivatives on Lower Circle
		else if (pId_d < xQuads)
		{
			deriv[0] = 0.0;
			deriv[1] = -dY;
			deriv[2] = 0.0;
		}

		// Derivatives on Straight-Up-Segment
		else if (pId_d < (xQuads + yQuads))
		{
			deriv[0] = 0.0;
			deriv[1] = 0.0;
			deriv[2] = dZ;
		}

		// Derivatives on Upper Circle
		else if (pId_d < (2 * xQuads + yQuads))
		{
			deriv[0] = 0.0;
			deriv[1] = dY;
			deriv[2] = 0.0;
		}

		// Derivatives on Straight-Down-Segment
		else if (pId_d < numPoints)
		{
			deriv[0] = 0.0;
			deriv[1] = 0.0;
			deriv[2] = -dZ;
		}


	derivatives->InsertNextTuple(deriv);
	}
	derivatives->Print(std::cout);

	inputPatch->GetPointData()->SetVectors(derivatives);

	// INSERT DERIVATIVES DONE.

	vtkDbiharStatic::WritePolyData(inputPatch, std::string(argv[0]) + "_patch.vtp");

	// vtkDbiharStatic::ShowPolyData(inputPatch); // Show Boundaries + Derivatives

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	// Set the bounds of the UV space.
	patchFilter->SetA(0.0);
	patchFilter->SetB(B);
	patchFilter->SetC(0.0);
	patchFilter->SetD(D);
	// Set the boundary conditions.
	patchFilter->SetMQuads(xQuads);
	patchFilter->SetNQuads(yQuads);
	// Set solution method.
	patchFilter->SetIFlag(2);

	patchFilter->SetInputData(inputPatch);

	patchFilter->Update();

	// patchFilter->Print(std::cout);

	vtkPolyData *outputPatch = patchFilter->GetOutput();

	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetDimensions(xQuads + 1, yQuads + 1, 1);
	structuredGrid->SetPoints(outputPatch->GetPoints());

	vtkDbiharStatic::ShowPolyDataWithGrid(inputPatch, structuredGrid); // Show Boundaries, Derivatives + Grid

	// Convert vtkStructuredGrid object to vtkPolydata.
	vtkSmartPointer<vtkStructuredGridGeometryFilter> gridGeometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	gridGeometryFilter->SetInputData(structuredGrid);
	gridGeometryFilter->Update();

	// vtkDbiharStatic::ShowPolyData(gridGeometryFilter->GetOutput()); // Show Grid
	vtkDbiharStatic::WritePolyData(gridGeometryFilter->GetOutput(), std::string(argv[0]) + "_solution.vtp");

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
