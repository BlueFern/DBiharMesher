/*
 * Patch filter for one branch (cylinder). Boundary condition on side of branch (outside)
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
	double L 		= 40;					//Set length
	double R 		= 10;					//Set radius
	int    uQuads 	= 24;					//Must be even! Number of quads along the length
	int    vQuads 	= 36;					//Must be even! Number of quads along the arch
	double angleBif	= 80;					//Angle of the bifurcation (in degrees)
	double deg2rad	= (vtkMath::Pi()/180);	//Degree to radiant converter

	double circle 	= 2*(vtkMath::Pi());	//Circle for the tube
	double derivL 	= 1.25*R;				//Set magnitude of derivatives along the length
	double derivA 	= L;					//Set magnitude of derivatives along the arch
	double alpha	= (angleBif/2)*deg2rad;			//Angle 1 for
	double beta		= ((360-angleBif)/4)*deg2rad;	//

	//define all further variables
	double point[3] = {0.0};
	double deriv[3] = {0.0};
	double gamma	= 0.0;
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
			x = (L-R/tan(beta))*((double)pId/uQuads);
			y = -R;
			z = 0;
			if(pId != 0)
			{
				dx = 0.0;
				dy = 0.0;
				dz = -derivL;
			}
		}
		//Set coordinates and derivatives for points on second segment (arch at bifurcation)
		else if(pId < (vQuads + uQuads))
		{
			gamma = circle*(((double)pId-uQuads)/vQuads);

			y = -R*cos(gamma);
			z = R*sin(gamma);

				if(gamma < (circle/4) || gamma > (3*circle/4))
				{
					x = L - (fabs(y)/tan(beta));
				}
				else
				{
					x = L - (fabs(y)/tan(alpha));
				}

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
			x = (L-R/tan(beta))*(1-((double)pId-(uQuads+vQuads))/uQuads);
			y = -R*cos(circle);
			z = R*sin(circle);

			if(pId != (vQuads + uQuads))
			{
				dx = 0.0;
				dy = 0.0;
				dz = derivL;
			}
		}
		//Set coordinates and derivatives for points on fourth segment (arch at x=0)
		else
		{
			gamma = circle*(1-(((double)pId-(2*uQuads+vQuads))/(double)vQuads));
			x = 0;
			y = -R*cos(gamma);
			z = R*sin(gamma);
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

	// vtkDbiharStatic::ShowPolyData(inputPatch);
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
