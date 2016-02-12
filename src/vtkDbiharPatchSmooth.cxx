#include <vtkAppendPolyData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkPolyLine.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkCellData.h>
#include <vtkObjectFactory.h>

#include "vtkDbiharPatchFilter.h"
#include "vtkDbiharPatchSmooth.h"
#include "vtkDbiharStatic.h"


vtkStandardNewMacro(vtkDbiharPatchSmooth);

vtkDbiharPatchSmooth::vtkDbiharPatchSmooth()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
	this->NumRadialQuads = 0;
}

int vtkDbiharPatchSmooth::RequestData(vtkInformation *vtkNotUsed(request),
		vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	// Test input data + params.
	if (this->NumRadialQuads == 0)
	{
		vtkErrorMacro("Must set the number of radial quads (halved)");
		exit(EXIT_FAILURE);
	}

	/*
	 * The following code only works under the following circumstances:
	 *
	 * -> there is only one bifurcation (ergo 3 branches)
	 * -> centreline is 2 dimensional in the x-y-plane
	 * -> parent branch grows in x-direction
	 * -> specific order of points must be given !?
	 */

	// Variables to get important points from the model
	int numPoints = input->GetNumberOfPoints() ;
	double *controlPoint;
	double xBif = 0.0;
	double endPoint1X = 0.0;		// x coordinate of end point daughter branch 1
	double endPoint1Y = 0.0;		// y coordinate of end point daughter branch 1
	double endPoint1Z = 0.0;		// z coordinate of end point daughter branch 1
	double endPoint2X = 0.0;		// x coordinate of end point daughter branch 2
	double endPoint2Y = 0.0;		// y coordinate of end point daughter branch 2
	double endPoint2Z = 0.0;		// z coordinate of end point daughter branch 2
	double xControl = 0.0;			// x coordinate of control point
	double yControl = 0.0;			// y coordinate of control point
	double zControl = 0.0;			// z coordinate of control point
	int xStep = 0;					// Size of step in x-direction, daughter branch 1
	int step = 0;					// Actual size of step in x-direction, daughter branch 1
	int c = 0;						// Control variable
	int pidRing = 0;				// Point ID (step = full ring)

	for (pidRing = 0; c < 2; pidRing = pidRing + 2*this->NumRadialQuads)
	{
		controlPoint = input->GetPoint(pidRing);
		xControl = controlPoint[0];
		yControl = controlPoint[1];
		zControl = controlPoint[2];

		if (yControl != 0.0 && c == 0)
		{
			xBif = input->GetPoint(pidRing - 2*this->NumRadialQuads)[0];
			endPoint2X = input->GetPoint(numPoints - 2*this->NumRadialQuads)[0];	// Get point coordinates for the last point of daughter branch 2
			endPoint2Y = input->GetPoint(numPoints - 2*this->NumRadialQuads)[1];
			endPoint2Z = input->GetPoint(numPoints - 2*this->NumRadialQuads)[2];
			xStep = round(xControl - xBif);
			c = 1;
		}
		if (c == 1)
		{
			controlPoint = input->GetPoint(pidRing + 2*this->NumRadialQuads);
			step = round(controlPoint[0] - xControl);

			if (step != xStep)
			{
				endPoint1X = xControl;	// Get point coordinates for the last point of daughter branch 1
				endPoint1Y = yControl;
				endPoint1Z = zControl;
				c = 2;
			}
		}
	}

	// Calculating characteristic properties
	double Length0 = xBif;													// Length of parent branch
	double Length1 = sqrt(pow(endPoint1X - xBif, 2) + pow(endPoint1Y, 2));	// Length of daughter branch 1
	double Length2 = sqrt(pow(endPoint2X-xBif, 2) + pow(endPoint2Y, 2));	// Length of daughter branch 2
	double Radius1 = sqrt(pow(endPoint1Z, 2));								// Radius of daughter branch 1 (calculated at end of daughter branch 1)
	double Radius2 = sqrt(pow(endPoint2Z, 2));								// Radius of daughter branch 2 (calculated at end of daughter branch 2)
	double Alpha = atan((endPoint1Y)/(endPoint1X-xBif));					// Calculating angle between X-axis and daughter branch 1
	double Beta = atan((-endPoint2Y)/(endPoint2X-xBif));					// Calculating angle between X-axis and daughter branch 2
	double Gamma = (Alpha + Beta)/2 - Beta;									// Calculating angle between X-axis and daughter-daughter-boundary

	// Calculating number of axial quads per branch
	vtkIdType pid = 0;
	int numAxialQuads0 = round(Length0/(vtkDbiharStatic::EC_AXIAL * 4));	// Number of ECs per quad length!
	int numAxialQuads1 = round(Length1/(vtkDbiharStatic::EC_AXIAL * 4));	// Number of ECs per quad length!
	int numAxialQuads2 = round(Length2/(vtkDbiharStatic::EC_AXIAL * 4));	// Number of ECs per quad length!

	if (numAxialQuads0 % 2) // Check if number of axial quads is even, subtract 1 if not (parent branch)
	{
		numAxialQuads0 = numAxialQuads0 - 1;
	}
	if (numAxialQuads1 % 2) // Check if number of axial quads is even, subtract 1 if not (daughter branch 1)
	{
		numAxialQuads1 = numAxialQuads1 - 1;
	}
	if (numAxialQuads2 % 2) // Check if number of axial quads is even, subtract 1 if not (daughter branch 2)
	{
		numAxialQuads2 = numAxialQuads2 - 1;
	}

	// Calculating number of points/cells per branch
	int numPointsBranch0 = 2*this->NumRadialQuads * numAxialQuads0 + 2*this->NumRadialQuads;	// Number of points parent branch
	int numPointsBranch1 = 2*this->NumRadialQuads * numAxialQuads1 + 2*this->NumRadialQuads;	// Number of points daughter branch 1
	int numPointsBranch2 = 2*this->NumRadialQuads * numAxialQuads2 + 2*this->NumRadialQuads;	// Number of points daughter branch 2
	int numCellsBranch0 = 2*this->NumRadialQuads * numAxialQuads0;								// Number of cells parent
	int numCellsBranch1 = 2*this->NumRadialQuads * numAxialQuads1;								// Number of cells daughter 1
	int numCellsBranch2 = 2*this->NumRadialQuads * numAxialQuads2;								// Number of cells daughter 2

	// Derivatives for the DbiharPatchFilter
	double dGROUND;					// Derivatives between daughter and parent
	double dBIF;					// Derivatives between daughter and daughter
	double dPEAK;					// Derivatives on peak (point of bifurcation)
	double dUP;						// Derivatives along daughter
	double dTOP;					// Derivatives on end of daughter

	double *exPoint;				// Declaration of 'extraction' point variable
	double deriv[3] = {0.0};		// Declaration of derivatives variable
	double Angle;					// Declaration of variable 'angle'

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();
	vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
	derivatives->SetName(vtkDbiharPatchFilter::DERIV_ARR_NAME);
	derivatives->SetNumberOfComponents(3);

	vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New(); // Adding the three branches to one model

	vtkSmartPointer<vtkPolyData> trunk = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> newCellArray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();

	// Keeping the parent branch
	for (int i = 0; i < numPointsBranch0; i++)
	{
		newPoints->InsertNextPoint(input->GetPoint(i));
	}
	for (int i = 0; i < numCellsBranch0; i++)
	{
		newCellArray->InsertNextCell(input->GetCell(i));
	}
	trunk->SetPoints(newPoints);
	trunk->SetPolys(newCellArray);
	appendPolyData->AddInputData(trunk);

	// Processing each daughter branch at a time
	for (int k = 0; k < 2; k++)
	{
		// Get points for boundary 1
		if(k == 0)
		{
			// Setting parameter for derivatives (daughter branch 1)
			dGROUND	= Length1;			// Derivatives between daughter1 and parent
			dBIF	= 0.25*Length1;		// Derivatives between daughter1 and daughter2
			dPEAK	= 1.70*Length1;		// Derivatives in point of bifurcation
			dUP		= 2*Radius1; 		// Derivatives on straight boundaries
			dTOP	= Length1;			// Derivatives on end of daughter branch
			Angle	= Alpha;			// Angle between parent and daughter1 axis

			// Get points for first boundary segment (circle at bifurcation)
			for(pid = numPointsBranch0 + this->NumRadialQuads/2; pid < numPointsBranch0 + 2*this->NumRadialQuads; pid++)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}
			for(pid = numPointsBranch0; pid < numPointsBranch0 + this->NumRadialQuads/2; pid++)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}
			// Get points for second boundary segment (straight line towards end of daughter branch)
			for(pid = numPointsBranch0 + this->NumRadialQuads/2; pid < numPointsBranch0 + numPointsBranch1 - 1.5*this->NumRadialQuads; pid = pid + 2*this->NumRadialQuads)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}
			//Get points for third boundary segment (circle at end of daughter branch)
			for(pid = numPointsBranch0 + numPointsBranch1 - 1.5*this->NumRadialQuads; pid >= numPointsBranch0 + numPointsBranch1 - 2*this->NumRadialQuads; --pid)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}
			for(pid = numPointsBranch0 + numPointsBranch1 - 1; pid > numPointsBranch0 + numPointsBranch1 - 1.5*this->NumRadialQuads; --pid)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}
			//Get points for fourth boundary segment (straight line towards parent branch)
			for(pid = numPointsBranch0 + numPointsBranch1 - 1.5*this->NumRadialQuads; pid > numPointsBranch0 + this->NumRadialQuads/2; pid = pid - 2*this->NumRadialQuads)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}

			// Setting derivatives for daughter branch 1
			for (int counter = 0; counter < points->GetNumberOfPoints(); counter++)
			{
				// Derivatives = 0 for corners
				if (counter == 0 || counter == 2*this->NumRadialQuads || counter == 2*this->NumRadialQuads + numAxialQuads1 || counter == 4*this->NumRadialQuads + numAxialQuads1)
				{
					deriv[0] = 0.0;
					deriv[1] = 0.0;
					deriv[2] = 0.0;
					derivatives->InsertNextTuple(deriv);
				}
				// Derivatives between daughter and parent branch (first half)
				else if (counter < this->NumRadialQuads/2)
				{
					deriv[0] = -cos(Angle * (0.5 + (double)counter/this->NumRadialQuads)) * dGROUND;
					deriv[1] = -sin(Angle * (0.5 + (double)counter/this->NumRadialQuads)) * dGROUND;
					deriv[2] = 0.0;
					derivatives->InsertNextTuple(deriv);
				}
				else if (counter > 1.5*this->NumRadialQuads && counter < 2*this->NumRadialQuads)
				{
					deriv[0] = -cos(Angle - (Angle*(counter-1.5*this->NumRadialQuads)/this->NumRadialQuads)) * dGROUND;
					deriv[1] = -sin(Angle - (Angle*(counter-1.5*this->NumRadialQuads)/this->NumRadialQuads)) * dGROUND;
					deriv[2] = 0.0;
					derivatives->InsertNextTuple(deriv);
				}
				// Derivatives for the peaks/tips
				else if (counter == this->NumRadialQuads/2 || counter == 1.5*this->NumRadialQuads)
				{
					deriv[0] = -cos(Angle) * dPEAK;
					deriv[1] = -sin(Angle) * dPEAK;
					deriv[2] = 0.0;
					derivatives->InsertNextTuple(deriv);
				}
				// Derivatives between daughter and daughter branch
				else if (counter < 1.5*this->NumRadialQuads)
				{
					if (counter <= this->NumRadialQuads)
					{
						deriv[0] =  sin(Gamma) * (dBIF*(2.0 - 2.0*(counter - (0.5*this->NumRadialQuads))/this->NumRadialQuads));
						deriv[1] = -cos(Gamma) * (dBIF*(2.0 - 2.0*(counter - (0.5*this->NumRadialQuads))/this->NumRadialQuads));
						deriv[2] = 0.0;
					}
					else
					{
						deriv[0] =  sin(Gamma) * (dBIF * (1.0 + 2.0*(counter - (this->NumRadialQuads))/this->NumRadialQuads));
						deriv[1] = -cos(Gamma) * (dBIF * (1.0 + 2.0*(counter - (this->NumRadialQuads))/this->NumRadialQuads));
						deriv[2] = 0.0;
					}
					derivatives->InsertNextTuple(deriv);
				}
				// Derivatives along circle (end of daughter branch)
				else if (counter > 2*this->NumRadialQuads + numAxialQuads1 && counter < 4*this->NumRadialQuads + numAxialQuads1)
				{
					deriv [0] = dTOP *cos(Angle);
					deriv [1] = dTOP *sin(Angle);
					deriv [2] = 0.0;
					derivatives->InsertNextTuple(deriv);
				}
				else
				// Derivatives along the branch (both directions)
				{
					deriv [0] = 0.0;
					deriv [1] = 0.0;
					deriv [2] = -dUP;
					if (counter < 4*this->NumRadialQuads + numAxialQuads1)
					{
						deriv [2] = -deriv [2];
					}
					derivatives->InsertNextTuple(deriv);
				}
			}
		}

		// Get points for boundary 2
		if(k == 1)
		{
			// Setting parameter for derivatives (daughter branch 2)
			dGROUND = Length2;			//Length 8800				// Derivatives between parent and daughter
			dBIF = 0.25*Length2;		//0.25*Length 2200			// Derivatives between daughter1 and daughter2
			dPEAK = 1.70*Length2;		//1.70*Length 14960			// Derivative in point of singularity
			dUP = 2*Radius2; 			//2*Radius 2200				// Derivatives on straight boundaries
			dTOP = Length2;				//Length 8800				// Derivatives on end of daughter branch
			Angle = Beta;				//							// Angle between parent and daughter2 axis

			// Get points for first boundary segment (circle at bifurcation)
			for(pid = numPointsBranch0 + numPointsBranch1 + 1.5*this->NumRadialQuads; pid >= numPointsBranch0 + numPointsBranch1; --pid)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}
			for(pid = numPointsBranch0 + numPointsBranch1 + 2*this->NumRadialQuads - 1; pid > numPointsBranch0 + numPointsBranch1 + 1.5*this->NumRadialQuads; --pid)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}
			// Get points for second boundary segment (straight line towards end of daughter branch)
			for(pid = numPointsBranch0 + numPointsBranch1 + 1.5*this->NumRadialQuads; pid < numPoints - this->NumRadialQuads/2; pid = pid + 2*this->NumRadialQuads)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}
			// Get points for third boundary segment (circle at end of daughter branch)
			for(pid = numPoints - this->NumRadialQuads/2; pid < numPoints; pid++)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}
			for(pid = numPoints - 2*this->NumRadialQuads; pid < numPoints - this->NumRadialQuads/2; pid++)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}
			// Get points for fourth boundary segment (straight line towards parent branch)
			for(pid = numPoints - this->NumRadialQuads/2; pid > numPointsBranch0 + numPointsBranch1 + 1.5*this->NumRadialQuads; pid = pid - 2*this->NumRadialQuads)
			{
				exPoint = input->GetPoint(pid);
				boundary->GetPointIds()->InsertNextId(points->InsertNextPoint(exPoint));
			}

			// Setting derivatives for daughter branch 2
			for (int counter = 0; counter < points->GetNumberOfPoints(); counter++)
			{
				// Derivatives = 0 for corners
				if (counter == 0 || counter == 2*this->NumRadialQuads || counter == 2*this->NumRadialQuads + numAxialQuads2 || counter == 4*this->NumRadialQuads + numAxialQuads2)
				{
					deriv[0] = 0.0;
					deriv[1] = 0.0;
					deriv[2] = 0.0;
					derivatives->InsertNextTuple(deriv);
				}
				// Derivatives between daughter and parent branch (first half)
				else if (counter < this->NumRadialQuads/2)
				{
					deriv[0] = -cos(Angle * (0.5 + (double)counter/this->NumRadialQuads)) * dGROUND;
					deriv[1] = sin(Angle * (0.5 + (double)counter/this->NumRadialQuads)) * dGROUND;
					deriv[2] = 0.0;
					derivatives->InsertNextTuple(deriv);
				}
				else if (counter > 1.5*this->NumRadialQuads && counter < 2*this->NumRadialQuads)
				{
					deriv[0] = -cos(Angle - (Angle*(counter-1.5*this->NumRadialQuads)/this->NumRadialQuads)) * dGROUND;
					deriv[1] = sin(Angle - (Angle*(counter-1.5*this->NumRadialQuads)/this->NumRadialQuads)) * dGROUND;
					deriv[2] = 0.0;
					derivatives->InsertNextTuple(deriv);
				}
				// Derivatives for the peaks/tips
				else if (counter == this->NumRadialQuads/2 || counter == 1.5*this->NumRadialQuads)
				{
					deriv[0] = -cos(Angle) * dPEAK;
					deriv[1] = sin(Angle) * dPEAK;
					deriv[2] = 0.0;
					derivatives->InsertNextTuple(deriv);
				}
				// Derivatives between daughter and daughter branch
				else if (counter < 1.5*this->NumRadialQuads)
				{
					if (counter <= this->NumRadialQuads)
					{
						deriv[0] = -sin(Gamma) * (dBIF*(2.0 - 2.0*(counter - (0.5*this->NumRadialQuads))/this->NumRadialQuads));
						deriv[1] =  cos(Gamma) * (dBIF*(2.0 - 2.0*(counter - (0.5*this->NumRadialQuads))/this->NumRadialQuads));
						deriv[2] = 0.0;
					}
					else
					{
						deriv[0] = -sin(Gamma) * (dBIF * (1.0 + 2.0*(counter - (this->NumRadialQuads))/this->NumRadialQuads));
						deriv[1] =  cos(Gamma) * (dBIF * (1.0 + 2.0*(counter - (this->NumRadialQuads))/this->NumRadialQuads));
						deriv[2] = 0.0;
					}
					derivatives->InsertNextTuple(deriv);
				}
				// Derivatives along circle (end of daughter branch)
				else if (counter > 2*this->NumRadialQuads + numAxialQuads2 && counter < 4*this->NumRadialQuads + numAxialQuads2)
				{
					deriv [0] = dTOP *cos(Angle);
					deriv [1] = -dTOP *sin(Angle);
					deriv [2] = 0.0;
					derivatives->InsertNextTuple(deriv);
				}
				else
				// Derivatives along the branch (both directions)
				{
					deriv [0] = 0.0;
					deriv [1] = 0.0;
					deriv [2] = -dUP;
					if (counter < 4*this->NumRadialQuads + numAxialQuads2)
					{
						deriv [2] = -deriv [2];
					}
					derivatives->InsertNextTuple(deriv);
				}
			}
		}

		vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
		boundaries->InsertNextCell(boundary);
		vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
		inputPatch->SetPoints(points);
		inputPatch->SetLines(boundaries);
		inputPatch->GetPointData()->SetVectors(derivatives);

		vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

		// Set the bounds of the UV space.
		patchFilter->SetA(0.0);
		patchFilter->SetB(3.0);
		patchFilter->SetC(0.0);
		patchFilter->SetD(1.0);
		// Set the boundary conditions.
		patchFilter->SetMQuads(2*this->NumRadialQuads);
		if (k == 0)
		{
			patchFilter->SetNQuads(numAxialQuads1);
		}
		else
		{
			patchFilter->SetNQuads(numAxialQuads2);
		}
		// Set solution method.
		patchFilter->SetIFlag(2);

		patchFilter->SetInputData(inputPatch);
		patchFilter->Update();

		//  patchFilter->Print(std::cout);

		vtkPolyData *outputPatch = patchFilter->GetOutput();

		vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
		if (k == 0)
		{
			structuredGrid->SetDimensions(2*this->NumRadialQuads + 1, numAxialQuads1 + 1, 1);
		}
		else
		{
			structuredGrid->SetDimensions(2*this->NumRadialQuads + 1, numAxialQuads2 + 1, 1);
		}
		structuredGrid->SetPoints(outputPatch->GetPoints());

//				 vtkDbiharStatic::ShowPolyDataWithGrid(inputPatch, structuredGrid);

		vtkSmartPointer<vtkStructuredGridGeometryFilter> gridGeometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
		gridGeometryFilter->SetInputData(structuredGrid);
		gridGeometryFilter->Update();

		appendPolyData->AddInputData(gridGeometryFilter->GetOutput());

#if 0
		// Write each branches boundary line and solution in separate files
		if(k == 0)
		{
			vtkDbiharStatic::WritePolyData(inputPatch, "boundary_1.vtp");
			vtkDbiharStatic::WritePolyData(gridGeometryFilter->GetOutput(), "solution_1.vtp");

		}
		else if(k == 1)
		{
			vtkDbiharStatic::WritePolyData(inputPatch, "boundary_2.vtp");
			vtkDbiharStatic::WritePolyData(gridGeometryFilter->GetOutput(), "solution_2.vtp");
		}
#endif
		inputPatch->Initialize();
		boundary->Initialize();
		points->Initialize();
		derivatives->Initialize();
	}

	appendPolyData->Update();

	vtkSmartPointer<vtkIntArray> branchIdCellData = vtkSmartPointer<vtkIntArray>::New();
	int branchId = 0;
	branchIdCellData->SetName(vtkDbiharStatic::BRANCH_ID_ARR_NAME);

	for (int cellId = 0; cellId < input->GetNumberOfCells(); cellId++)
	{
		if (cellId >= numCellsBranch0 && cellId < numCellsBranch0 + numCellsBranch1)
		{
			branchId = 1;
		}
		if (cellId >= numCellsBranch0 + numCellsBranch1)
		{
			branchId = 2;
		}
		branchIdCellData->InsertNextValue(branchId);
	}

	appendPolyData->GetOutput()->GetCellData()->SetScalars(branchIdCellData);
	output->DeepCopy(appendPolyData->GetOutput());
	return 1;
}
void vtkDbiharPatchSmooth::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Number of radial quads: " << this->NumRadialQuads << "\n";
}
