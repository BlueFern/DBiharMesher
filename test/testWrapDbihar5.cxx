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
#include "showPolyData.h"

#include "wrapDbiharConfig.h"

double radToDeg(double angleInRad)
{
	return angleInRad * (180.0/vtkMath::Pi());
}

double quadraticFunction(double x)
{
	// A and B calculator from three points: http://www.softschools.com/math/algebra/quadratic_functions/quadratic_function_with_three_points/
	// Points [0,1], [0.5,0.8], [1,1].
	double a = 3.1999;
	double b = -1.5999;
	double c = 1;
	return a * x * x + b * x + c;
}

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	// Create the sides of a patch as 3D points.
	double alpha = vtkMath::Pi() / 4.0;
	double x = 0.0;
	double y = 0.0;
	double z1 = -10.0;
	double z2 = 10.0;
	double halfLen = 60.0;
	double arc = vtkMath::Pi();

	int cQuads = 18; // m = 17. Num quads should be even, to make sure m is odd.
	int yQuads = 60; // n = 59. Num quads should be even, to make sure n is odd.

	vtkIdType pIds = (cQuads + yQuads) * 2;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();

	// Insert boundary points. The boundary has four segments.
	// The coordinates of the points are calculated specific to the current boundary segment.
	double radius = (z2 - z1) / 2;
	for(vtkIdType pId = 0; pId < pIds; pId++)
	{
		// This point is declared inside the loop to make sure it is (0,0,0) at the start of every iteration.
		double point[3] = {0.0};

		// Terminal arc.
		// Inserting points along the y = y1 boundary segment.
		if(pId < cQuads)
		{
			// First, rotate radius vector around (0,0,0).
			double dA = pId / (double)cQuads;
			double tmpPoint[3] = {0.0};
			if(dA < 0.5)
			{
				double angle = arc * dA;
				tmpPoint[1] = sin(angle) * radius;
				tmpPoint[2] = -cos(angle) * radius;
			}
			else
			{
				double angle = arc * dA - vtkMath::Pi() / 2.0;
				tmpPoint[1] = cos(angle) * radius;
				tmpPoint[2] = sin(angle) * radius;
			}
			//PrintPoint(tmpPoint); std::cout << std::endl;
			// Then rotate it into the plane of the terminal arc, and translate to the end/start point of this side.
			point[0] = -cos(alpha) * halfLen + cos(alpha) * tmpPoint[1];
			point[1] = sin(alpha) * halfLen + sin(alpha) * tmpPoint[1];
			point[2] = tmpPoint[2];
		}
		// Line parallel to the centreline.
		// Inserting points along the x = x2 boundary segment.
		else if(pId < cQuads + yQuads)
		{
			// Parametric coordinate along this edge: [0,1].
			double dL = (pId - cQuads) / (double)yQuads;
			if(dL < 0.5)
			{
				point[0] = x - (halfLen * 2.0 * std::fabs(dL - 0.5)) * cos(alpha);
				point[1] = y + (halfLen * 2.0 * std::fabs(dL - 0.5)) * sin(alpha);
			}
			else
			{
				point[0] = x + (halfLen * 2.0 * (dL - 0.5)) * cos(alpha);
				point[1] = y + (halfLen * 2.0 * (dL - 0.5)) * sin(alpha);
			}
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
			point[0] = cos(alpha) * halfLen - cos(alpha) * tmpPoint[1];
			point[1] = sin(alpha) * halfLen + sin(alpha) * tmpPoint[1];
			point[2] = tmpPoint[2];
		}
		// Line parallel to the centreline.
		// Inserting points along the x = x1 boundary segment.
		else
		{
			// Parametric coordinate along this edge: [0,1].
			double dL = (pId - cQuads - yQuads - cQuads) / (double)yQuads;
			if(dL < 0.5)
			{
				point[0] = x + (halfLen * 2.0 * std::fabs(dL - 0.5)) * cos(alpha);
				point[1] = y + (halfLen * 2.0 * std::fabs(dL - 0.5)) * sin(alpha);
			}
			else
			{
				point[0] = x - (halfLen * 2.0 * (dL - 0.5)) * cos(alpha);
				point[1] = y + (halfLen * 2.0 * (dL - 0.5)) * sin(alpha);
			}
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

	// Base magnitude of derivatives on axial boundaries.
	double magAx = 60.0; // Should be calculated from radius.
	// Base magnitude of derivatives on terminal boundaries.
	double magTerm = 35; // Should be calculated from radius.
	// Additional scaling for derivatives on axial boundaries.
	double dMag = 0.8;
	double dAlpha = 0.15;
	// Rotational scaling for derivatives.
	//double dAlpha = 1.5;
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
			// Skipping first point.
			if(pId != 0)
			{
				deriv[0] = -cos(alpha) * magTerm;
				deriv[1] = sin(alpha) * magTerm;
			}
		}
		// Inserting derivatives along the x = x2 boundary segment, skipping the corner case.
		else if(pId < cQuads + yQuads)
		{
			// Skipping first point.
			if(pId != cQuads)
			{
				// Parametric position along this edge: (0, 1).
				double dL = (pId - cQuads) / double(yQuads);
				if(dL < 0.5)
				{
					// Close to 0 the angle of the derivative is nearly pi - alpha; close to 0.5 the angle is -2 * alpha.
					// It is simpler to calculate the derivative from alpha to pi / 2 and the flip the vector.
					double derivAlpha = alpha + (vtkMath::Pi() - 2.0 * alpha) * dL + dAlpha * 2.0 * dL;
					// Scaling is from 1 to 1 + dMag.
					deriv[0] = (-cos(derivAlpha) * magAx) * (1.0 + dMag * dL * 2.0) * quadraticFunction(dL);
					deriv[1] = (-sin(derivAlpha) * magAx) * (1.0 + dMag * dL * 2.0) * quadraticFunction(dL);
				}
				else
				{
					// Close to 0.5 the angle of the derivative is nearly -2 * alpha; close to 1 the angle is -alpha.
					// It is simpler to calculate the derivative from pi / 2 to alpha and then flip the y component of the vector.
					double derivAlpha = vtkMath::Pi() / 2.0 - (alpha * 2.0 * (dL - 0.5));
					if(dL > 0.5)
					{
						derivAlpha += dAlpha * 2.0 * std::fabs(dL - 1.0);
					}
					// Scaling is from 1 + dMag to 1.
					deriv[0] = (cos(derivAlpha) * magAx) * (1.0 + dMag * std::fabs(dL - 1.0) * 2.0) * quadraticFunction(dL - 0.5);
					deriv[1] = (-sin(derivAlpha) * magAx) * (1.0 + dMag * std::fabs(dL - 1.0) * 2.0) * quadraticFunction(dL - 0.5);
				}
			}
		}
		// Inserting derivatives along the y = y2 boundary segment, skipping the corner case.
		else if(pId < cQuads * 2 + yQuads)
		{
			if(pId != cQuads + yQuads)
			{
				deriv[0] = cos(alpha) * magTerm;
				deriv[1] = sin(alpha) * magTerm;
			}
		}
		// Inserting derivatives along the x = x1 boundary segment, skipping the corner case.
		else
		{
			// Skipping the first point.
			if(pId != cQuads * 2 + yQuads)
			{
				// Parametric position along this edge: (0, 1).
				double dL = (pId - (cQuads * 2 + yQuads)) / double(yQuads);
				if(dL < 0.5)
				{
					// Close to 0 the angle of the derivative is nearly -alpha; close to 0.5 the angle is -pi / 2.
					// It is simpler to calculate the derivative from alpha to pi / 2 and flip the the x component of the vector.
					double derivAlpha = alpha + (vtkMath::Pi() - 2.0 * alpha) * dL + dAlpha * 2.0 * dL;
					// Scaling is from 1 to 1 + dMag.
					deriv[0] = (cos(derivAlpha) * magAx) * (1.0 + dMag * dL * 2.0) * quadraticFunction(dL);
					deriv[1] = (-sin(derivAlpha) * magAx) * (1.0 + dMag * dL * 2.0) * quadraticFunction(dL);
				}
				else
				{
					// Close to 0.5 the angle of derivative is nearly -pi / 2; close to 1 the angle is pi - alpha.
					// It is simpler to calculate the derivative from pi / 2 to alpha and then flip the the vector.
					double derivAlpha = vtkMath::Pi() / 2.0 - (alpha * 2.0 * (dL - 0.5));
					if(dL > 0.5)
					{
						derivAlpha += dAlpha * 2.0 * std::fabs(dL - 1.0);
					}
					// Scaling is from 1 + dMag to 1.
					deriv[0] = (-cos(derivAlpha) * magAx) * (1.0 + dMag * std::fabs(dL - 1.0) * 2.0) * quadraticFunction(dL - 0.5);
					deriv[1] = (-sin(derivAlpha) * magAx) * (1.0 + dMag * std::fabs(dL - 1.0) * 2.0) * quadraticFunction(dL - 0.5);
				}
			}
		}
		derivatives->InsertNextTuple(deriv);
	}

	inputPatch->GetPointData()->SetVectors(derivatives);

	showPolyData(inputPatch, NULL);

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

	showPolyData(inputPatch, structuredGrid);

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
