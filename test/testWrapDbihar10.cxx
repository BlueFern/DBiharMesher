/*
 * Patch filter for one branch (cylinder). Boundary condition on top of branch (upside)
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

/*
	//Defining the parameters
	double L 		= 40;					//Set length
	double R 		= 10;					//Set radius
	int    uQuads 	= 24;					//Must be even! Number of quads along the length
	int    vQuads 	= 36;					//Must be even! Number of quads along the arch
	double angleBif	= 80;					//Angle of the bifurcation (in degrees)
	double deg2rad	= (vtkMath::Pi()/180);	//

	double circle 	= 2*(vtkMath::Pi());	//Set angle of arch
	double derivL 	= 1.25*R;				//Set magnitude of derivatives along the length
	double derivA 	= L;					//Set magnitude of derivatives along the arch
	double alpha	= (angleBif/2)*deg2rad;
	double beta		= ((360-angleBif)/4)*deg2rad;

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
		//Set coordinates and derivatives for points on second segment (arch at bifurcation)
		else if(pId < (vQuads + uQuads))
		{
			gamma = circle*(((double)pId-uQuads)/vQuads);

			y = -R*cos(gamma);
			z = R*sin(gamma);

				if(gamma <= (circle/2))
				{
					x = L - (z/tan(alpha));
				}
				else if(gamma > (circle/2))
				{
					x = L - (fabs(z)/tan(beta));
				}

			if(pId != uQuads)
			{
				dx = 0.5*derivA;
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
*/


	double Length = 8800;
	double Radius = 1273.24;
	double dA = 30000;
	double dL = 2500;
	int circQuads = 20;
	int axiQuads = 34;
	int numQuads = 2*(circQuads + 2*axiQuads);
	double point [3] = {0.0};
	double deriv [3] = {0.0};
	double Angle = 30*(vtkMath::Pi()/180);
	double Step = 0.0;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();
	vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
	derivatives->SetName(vtkDbiharPatchFilter::DERIV_ARR_NAME);
	derivatives->SetNumberOfComponents(3);

	for (vtkIdType pId = 0; pId < numQuads; pId++)
	{
		deriv[0] = 0.0;
		deriv[1] = 0.0;
		deriv[2] = 0.0;

		if (pId < circQuads)
		{
			Step = (double)pId/circQuads;
			point[0] = -Length;
			point[1] = sin(vtkMath::Pi()*Step) * Radius;
			point[2] = cos(vtkMath::Pi()*Step) * Radius;

			if (pId != 0)
			{
				deriv[0] = -dA;
				deriv[1] = 0.0;
				deriv[2] = 0.0;
			}
		}

		else if (pId < (circQuads + axiQuads))
		{
			Step = ((double)pId-circQuads)/axiQuads;
			point[0] = Length * (Step - 1);
			point[1] = 0.0;
			point[2] = - Radius;

			if (pId != circQuads)
			{
				deriv[0] = 0.0;
				deriv[1] = -dL;
				deriv[2] = 0.0;
			}
		}

		else if (pId < circQuads + 2*axiQuads)
		{
			Step = ((double)pId-(circQuads+axiQuads))/axiQuads;
			point[0] = cos(Angle) * (Length * Step);
			point[1] = sin(Angle) * (Length * Step);
			point[2] = - Radius;

			if (pId == (circQuads + axiQuads))
			{
				deriv[0] = sin(Angle/2)*dL;
				deriv[1] = -cos(Angle/2)*dL;
				deriv[2] = 0.0;
			}
			else
			{
				deriv[0] = sin(Angle)*dL;
				deriv[1] = -cos(Angle)*dL;
				deriv[2] = 0.0;
			}
		}

		else if (pId < 2*(circQuads + axiQuads))
		{
			Step = ((double)pId-(circQuads + 2*axiQuads))/circQuads;
			point[0] = cos(Angle) * Length - (Radius * sin(vtkMath::Pi() * Step))*sin(Angle);
			point[1] = (Length * sin(Angle)) + (Radius * sin(vtkMath::Pi() * Step)*cos(Angle));
			point[2] = -cos(vtkMath::Pi() * Step) * Radius;

			if (pId != (circQuads + 2*axiQuads))
			{
				deriv[0] = cos(Angle)*dA;
				deriv[1] = sin(Angle)*dA;
				deriv[2] = 0.0;
			}
		}

		else if (pId < (2*circQuads + 3*axiQuads))
		{
			Step = ((double)pId - (2*circQuads + 2*axiQuads))/axiQuads;
			point[0] = cos(Angle) * Length - cos(Angle) * Length * Step;
			point[1] = sin(Angle) * Length - sin(Angle) * Length * Step;
			point[2] = Radius;

			if (pId != 2*(circQuads + axiQuads))
			{
				deriv[0] = sin(Angle)*dL;
				deriv[1] = -cos(Angle)*dL;
				deriv[2] = 0.0;
			}
		}

		else
		{
			Step = ((double)pId - (2*circQuads + 3*axiQuads))/axiQuads;
			point[0] = - Length * Step;
			point[1] = 0.0;
			point[2] = Radius;

			if (pId == (2*circQuads + 3*axiQuads))
			{
				deriv[0] = sin(Angle/2)*dL;
				deriv[1] = -cos(Angle/2)*dL;
				deriv[2] = 0.0;
			}
			else
			{
				deriv[0] = 0.0;
				deriv[1] = -dL;
				deriv[2] = 0.0;
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


//	vtkDbiharStatic::ShowPolyData(inputPatch);
	vtkDbiharStatic::WritePolyData(inputPatch, "sattle_test_inputPatch.vtp");

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	// Set the bounds of the UV space.
	patchFilter->SetA(0.0);
	//patchFilter->SetB(2.0/3.0);
	patchFilter->SetB(2.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(1.0);
	// Set the boundary conditions.
	patchFilter->SetMQuads(circQuads);
	patchFilter->SetNQuads(2*axiQuads);
	// Set solution method.
	patchFilter->SetIFlag(2);

	patchFilter->SetInputData(inputPatch);
	patchFilter->Update();

	//  patchFilter->Print(std::cout);

	vtkPolyData *outputPatch = patchFilter->GetOutput();

	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetDimensions(circQuads + 1, 2*axiQuads + 1, 1);
	structuredGrid->SetPoints(outputPatch->GetPoints());

	vtkDbiharStatic::ShowPolyDataWithGrid(inputPatch, structuredGrid);

	vtkSmartPointer<vtkStructuredGridGeometryFilter> gridGeometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	gridGeometryFilter->SetInputData(structuredGrid);
	gridGeometryFilter->Update();

	vtkDbiharStatic::WritePolyData(gridGeometryFilter->GetOutput(), "sattle_test_solution.vtp");


	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
