/*
 * Patch filter for cylinder fractions
 */

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

	//Defining the parameters
	double frac		= .75;					//Fraction of cylinder [0...1]
	double L 		= 8800;					//Set length
	double R 		= 1273.24;					//Set radius
	int    uQuads 	= 40;					//Must be even! Number of quads along the length
	int    vQuads 	= 34;					//Must be even! Number of quads along the arch

	double circle 	= frac*2*(vtkMath::Pi());	//Set angle of arch
	double derivL 	= (frac+0.25)*R;		//Set magnitude of derivatives along the length
	double derivA 	= L;					//Set magnitude of derivatives along the arch

	//define all further variables
	double point[3] = {0.0};
	double deriv[3] = {0.0};
	double alpha	= 0.0;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
	double dx = 0.0;
	double dy = 0.0;
	double dz = 0.0;
	vtkIdType pIds = (uQuads + vQuads) * 2;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();
	vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
	derivatives->SetName(vtkDbiharPatchFilter::DERIV_ARR_NAME);
	derivatives->SetNumberOfComponents(3);

	for(vtkIdType pId = 0; pId < pIds; pId++)
	{
		//Null derivatives for the corner cases
		dx = 0.0;
		dy = 0.0;
		dz = 0.0;
		//Set coordinates and derivatives for points on first segment (straight line at y=-R)
		//Starting at x=0, y=-R, z=0)
		if(pId < uQuads)
		{
			x = L*((double)pId/uQuads);
			y = -R;
			z = 0;
			if(pId != 0)
			{
				dx = 0.0;
				dy = 0.0;
				dz = -derivL;
			}
		}
		//Set coordinates and derivatives for points on second segment (arch at x=L)
		else if(pId < (vQuads + uQuads))
		{
			alpha = circle*(((double)pId-uQuads)/vQuads);
			x = L;
			y = -R*cos(alpha);
			z = R*sin(alpha);
			if(pId != uQuads)
			{
				dx = derivA;
				dy = 0.0;
				dz = 0.0;
			}
		}
		//Set coordinates and derivatives for points on third segment (straight line)
		else if(pId < (2*uQuads + vQuads))
		{
			x = L*(1-((double)pId-(uQuads+vQuads))/uQuads);
			y = -R*cos(circle);
			z = R*sin(circle);
			if(pId != (vQuads + uQuads))
			{
				dx = 0.0;
				dy = derivL*sin(circle);
				dz = derivL*cos(circle);
			}
		}
		//Set coordinates and derivatives for points on fourth segment (arch at x=0)
		else
		{
			alpha = circle*(1-(((double)pId-(2*uQuads+vQuads))/(double)vQuads));
			x = 0;
			y = -R*cos(alpha);
			z = R*sin(alpha);
			if(pId != (2*uQuads + vQuads))
			{
				dx = -derivA;
				dy = 0.0;
				dz = 0.0;
			}
		}
		point[0] = x;
		point[1] = y;
		point[2] = z;

		deriv[0] = dx;
		deriv[1] = dy;
		deriv[2] = dz;

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
	//patchFilter->SetB(2.0/3.0);
	patchFilter->SetB(1.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(5.0);
	// Set the boundary conditions.
	patchFilter->SetMQuads(uQuads);
	patchFilter->SetNQuads(vQuads);
	// Set solution method.
	patchFilter->SetIFlag(2);

	patchFilter->SetInputData(inputPatch);
	patchFilter->Update();

	//  patchFilter->Print(std::cout);

	vtkPolyData *outputPatch = patchFilter->GetOutput();

	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetDimensions(uQuads + 1, vQuads + 1, 1);
	structuredGrid->SetPoints(outputPatch->GetPoints());

	vtkDbiharStatic::ShowPolyDataWithGrid(inputPatch, structuredGrid);

	vtkSmartPointer<vtkStructuredGridGeometryFilter> gridGeometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	gridGeometryFilter->SetInputData(structuredGrid);
	gridGeometryFilter->Update();

	vtkDbiharStatic::WritePolyData(gridGeometryFilter->GetOutput(), std::string(argv[0]) + "_solution.vtp");


	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
