/*
 * 1/6th of a Bifurcation Patch Generation: Half of the Saddle.
 */

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

double quadraticFunction(double x)
{
	// A and B calculator from three points: http://www.softschools.com/math/algebra/quadratic_functions/quadratic_function_with_three_points/
	// Points [0,1], [0.5,0.8], [1,1].
	double a = 0.7999;
	double b = -0.7999;
	return a * x * x + b * x + 1;
}

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	// Create the sides of a patch as 3D points.
	double alpha = vtkMath::Pi() / 4.0;
	double x = 0.0;
	double y = 0.0;
	double z1 = -10.0;
	double z2 = 10.0;
	double Len = 60.0;
	double arc = vtkMath::Pi();

	int cQuads = 18; // m = 17. Num quads should be even, to make sure m is odd.
	int yQuads = 30; // n = 29. Num quads should be even, to make sure n is odd.

	vtkIdType pIds = (cQuads + yQuads) * 2;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();

	// Insert boundary points. The boundary has four segments.
	// The coordinates of the points are calculated specific to the current boundary segment.
	double radius = (z2 - z1) / 2;
	double dXNeckScale = 1.8;
	for(vtkIdType pId = 0; pId < pIds; pId++)
	{
		// This point is declared inside the loop to make sure it is (0,0,0) at the start of every iteration.
		double point[3] = {0.0};

		// Bifurcation arc.
		// Inserting points along the y = y1 boundary segment.
		// Scaling in the Y dimension around this edge is done for smoothing the surface when two patches are joined together.
		if(pId < cQuads)
		{
			double dA = pId / (double)cQuads;
			if(dA < 0.5)
			{
				double angle = arc * dA;
				point[1] = (sin(angle) * radius) * dXNeckScale;
				point[2] = -cos(angle) * radius;
			}
			else
			{
				double angle = arc * dA - vtkMath::Pi() / 2.0;
				point[1] = (cos(angle) * radius) * dXNeckScale;
				point[2] = sin(angle) * radius;
			}
			point[0] = x;
		}
		// Line parallel to the centreline.
		// Inserting points along the x = x2 boundary segment.
		else if(pId < cQuads + yQuads)
		{
			double l = (pId - cQuads) * (Len / (double)yQuads);
			point[0] = x + l * cos(alpha);
			point[1] = y + l * sin(alpha);
			point[2] = z2;
		}
		// Terminal arc.
		// Inserting points along the y = y2 boundary segment.
		else if(pId < cQuads * 2 + yQuads)
		{
			// First, rotate radius vector around (0,0,0).
			double dA = (pId - cQuads - yQuads) / (double)cQuads;
			double tmpPoint[3] = {0.0};
			if(dA < 0.5)
			{
				double angle = arc * dA;
				tmpPoint[1] = sin(angle) * radius;
				tmpPoint[2] = cos(angle) * radius;
			}
			else
			{
				double angle = arc * dA - vtkMath::Pi() / 2.0;
				tmpPoint[1] = cos(angle) * radius;
				tmpPoint[2] = -sin(angle) * radius;
			}
			// Then rotate it into the plane of the terminal arc, and translate to the end/start point of this side.
			point[0] = cos(alpha) * Len - cos(alpha) * tmpPoint[1];
			point[1] = sin(alpha) * Len + sin(alpha) * tmpPoint[1];
			point[2] = tmpPoint[2];
		}
		// Line parallel to the centreline.
		// Inserting points along the x = x1 boundary segment.
		else
		{
			double l = (pId - cQuads - yQuads - cQuads) * (Len / (double)yQuads);
			point[0] = cos(alpha) * Len - l * cos(alpha);
			point[1] = sin(alpha) * Len - l * sin(alpha);
			point[2] = z1;
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

	// Derivatives.
	vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
	derivatives->SetName(vtkDbiharPatchFilter::DERIV_ARR_NAME);
	derivatives->SetNumberOfComponents(3);

	// Base magnitude of derivatives on straight edges.
	double mag = 60.0; // Should be calculated from radius.
	// Additional scaling for derivatives on straight edges.
	double dMag = 0.70; // Should be calculated from magnitude.
	// Rotational scaling for derivatives.
	double dAlpha = 1.2;
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
				deriv[0] = -cos(alpha) * 10.0;
				deriv[1] = -sin(alpha) * 10.0;
				deriv[2] = 0.0;
			}
		}
		// Inserting derivatives along the x = x2 boundary segment, skipping the corner case.
		else if(pId < cQuads + yQuads)
		{
			// Parametric position along this edge: (0, 1).
			double dL = (pId - cQuads) / double(yQuads);
			// Close to 0 the angle of derivative is nearly alpha; close to 1 the angle of derivative is -alpha.
			double derivAlpha = alpha - (-dAlpha * dL + dAlpha) * alpha;
			if(pId != cQuads)
			{
				deriv[0] = (sin(derivAlpha) * mag) * (1 + (-dMag * dL + dMag)) * quadraticFunction(dL);
				deriv[1] = (-cos(derivAlpha) * mag) * (1 + (-dMag * dL + dMag)) * quadraticFunction(dL);
				deriv[2] = 0.0;
			}
		}
		// Inserting derivatives along the y = y2 boundary segment, skipping the corner case.
		else if(pId < cQuads * 2 + yQuads)
		{
			if(pId != cQuads + yQuads)
			{
				deriv[0] = cos(alpha) * 15.0;
				deriv[1] = sin(alpha) * 15.0;
				deriv[2] = 0.0;
			}
		}
		// Inserting derivatives along the x = x1 boundary segment, skipping the corner case.
		else
		{
			// Parametric position along this edge: (0, 1).
			double dL = (pId - (cQuads * 2 + yQuads)) / double(yQuads);
			// Close to 0 the angle of derivative is nearly alpha; close to 1 the angle of derivative is -alpha.
			double derivAlpha = alpha + (-dAlpha * dL) * alpha;
			if(pId != cQuads * 2 + yQuads)
			{
				deriv[0] = (sin(derivAlpha) * mag) * (1 + (dMag * dL)) * quadraticFunction(dL);
				deriv[1] = (-cos(derivAlpha) * mag) * (1 + (dMag * dL)) * quadraticFunction(dL);
				deriv[2] = 0.0;
			}
		}
		derivatives->InsertNextTuple(deriv);
	}

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
