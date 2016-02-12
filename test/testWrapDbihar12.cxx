/*
 * Patch filter for one branch (bifurcation)
 * First daughter branch
 * For full bifurcation 'testWrapDbihar13' and 'testWrapDbihar14' are needed
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
	double L 		= 8800;					//Set length
	double R 		= 1100;		// 1273.24	//Set radius
	int    uQuads 	= 40;					//Must be even! Number of quads along the length (axial)
	int    vQuads 	= 34;					//Must be even! Number of quads along the arch (circumferential)
	double angleBif	= 60;					//Angle of the bifurcation (in degrees)
	double depth1	= 1100;	// 1430.45		//Depth in Z direction
	double depth2	= 1485.01;	// 1485.01		//Depth in -Z direction
	double derivB	= 0.25*L;				//Set magnitude of derivatives along fist half of bifurcation
	double derivC	= 0.5*L;				//Set magnitude of derivatives along second half of bifurcation

	double deg2rad	= (vtkMath::Pi()/180);	//
	double circle 	= 2*(vtkMath::Pi());	//Set angle of arch
	double derivL 	= 1.25*R;				//Set magnitude of derivatives along the length
	double derivA 	= L;					//Set magnitude of derivatives along the circle (arch)
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
		//Set coordinates and derivatives for points on first segment (circle at x=-L)
		if(pId < uQuads)
		{
			gamma = circle*((double)pId/(double)uQuads);
			x = -L;
			y = -R*cos(gamma);
			z = -R*sin(gamma);
			if(pId != 0)
			{
				dx = -derivA;
				dy = 0.0;
				dz = 0.0;
			}
		}
		//Set coordinates and derivatives for points on second segment (straight line at y=-R)
		else if(pId < (uQuads + vQuads))
		{
			x = L*((((double)pId-uQuads)/vQuads)-1);
			y = -R;
			z = 0;
			if(pId != uQuads)
			{
				dx = 0.0;
				dy = 0.0;
				dz = -derivL;
			}
		}

		//Set coordinates and derivatives for points on third segment (two arches at bifurcation)
		//y and z calculated from equation for ellipse: (y^2)/(R^2) + (z^2)/(depth^2) = 1
		else if(pId < (2*uQuads + vQuads))
		{
					gamma = circle*(((double)pId- (uQuads+vQuads)) / uQuads);

						if(gamma <= (circle/4))
						{
							y = - (depth1 * R			  )/sqrt(pow(depth1,2) + pow(R,2)*pow(tan(gamma),2));
							z = + (depth1 * R * tan(gamma))/sqrt(pow(depth1,2) + pow(R,2)*pow(tan(gamma),2));
							x =  - (z/tan(alpha));

							if(pId != uQuads + vQuads)
							{
								dx = derivB*sin(alpha);
								dy = 0.0;
								dz = derivB*cos(alpha);
							}
						}
						else if(gamma <= (circle/2))
						{
							y = + (depth1 * R			  )/sqrt(pow(depth1,2) + pow(R,2)*pow(tan(gamma),2));
							z = - (depth1 * R * tan(gamma))/sqrt(pow(depth1,2) + pow(R,2)*pow(tan(gamma),2));
							x = - (z/tan(alpha));

							if(gamma == circle/2)
							{
								dx = derivC + derivB;
								dy = 0.0;
								dz = 0.0;
							}
							else
							{
								dx = derivB*sin(alpha);
								dy = 0.0;
								dz = derivB*cos(alpha);
							}
						}
						else if(gamma <= (circle*3/4))
						{
							y = + (depth2 * R			  )/sqrt(pow(depth2,2) + pow(R,2)*pow(tan(gamma),2));
							z = - (depth2 * R * tan(gamma))/sqrt(pow(depth2,2) + pow(R,2)*pow(tan(gamma),2));
							x = - (fabs(z)/tan(beta));

							dx = derivC*sin(beta);
							dy = 0.0;
							dz = -derivC*cos(beta);
						}
						else
						{
							y = - (depth2 * R			  )/sqrt(pow(depth2,2) + pow(R,2)*pow(tan(gamma),2));
							z = + (depth2 * R * tan(gamma))/sqrt(pow(depth2,2) + pow(R,2)*pow(tan(gamma),2));
							x = - (fabs(z)/tan(beta));

							dx = derivC*sin(beta);
							dy = 0.0;
							dz = -derivC*cos(beta);
						}
				}

		//Set coordinates and derivatives for points on fourth segment (straight line at y=-R)
		else
		{
					x = -L*((((double)pId-(2*uQuads+vQuads))/vQuads));
					y = -R;
					z = 0.0;

					if(pId != (2*uQuads + vQuads))
					{
						dx = 0.0;
						dy = 0.0;
						dz = derivL;
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
	patchFilter->SetB(3.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(1.0);
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
