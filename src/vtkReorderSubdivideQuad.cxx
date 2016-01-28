#include <algorithm>

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkCallbackCommand.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>

#include "vtkSubdivideQuadBrick.h"
#include "vtkReorderSubdivideQuad.h"

vtkStandardNewMacro(vtkReorderSubdivideQuad);

vtkReorderSubdivideQuad::vtkReorderSubdivideQuad()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Columns = 0;
	this->Rows = 0;
	this->Rotations = 0;
}

int vtkReorderSubdivideQuad::RequestData(vtkInformation *vtkNotUsed(request),
		vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	// Reorder the points (to rotate the quad).

	vtkSmartPointer<vtkIdList> pointsList = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	pointsList->SetNumberOfIds(4);
	points->SetNumberOfPoints(4);
	points->SetPoint(0,input->GetPoint((0 + (4 - this->Rotations)) % 4));
	points->SetPoint(1,input->GetPoint((1 + (4 - this->Rotations)) % 4));
	points->SetPoint(2,input->GetPoint((2 + (4 - this->Rotations)) % 4));
	points->SetPoint(3,input->GetPoint((3 + (4 - this->Rotations)) % 4));

	vtkSmartPointer<vtkPolyData> pointsData = vtkSmartPointer<vtkPolyData>::New();
	pointsData->SetPoints(points);

	vtkSmartPointer<vtkSubdivideQuadBrick> subdivideQuadBrick = vtkSmartPointer<vtkSubdivideQuadBrick>::New();
	subdivideQuadBrick->SetInputData(pointsData);


	//TODO: doc
	if (this->Rotations % 2 == 0) // EC
	{
		subdivideQuadBrick->SetColumns(this->Columns);
		subdivideQuadBrick->SetRows(this->Rows);
	}
	else // SMC
	{
		subdivideQuadBrick->SetColumns(this->Rows);
		subdivideQuadBrick->SetRows(this->Columns);

	}
	subdivideQuadBrick->Update();

	vtkSmartPointer<vtkPolyData> subdividedQuad = subdivideQuadBrick->GetOutput();
	subdividedQuad->BuildCells();

	vtkSmartPointer<vtkCellArray> cells = subdividedQuad->GetPolys();
	vtkSmartPointer<vtkPoints> oldPoints = subdividedQuad->GetPoints();

	// Now reorder cells in the subdivided quad - again necessary because it was initally rotated.
	// Note: Columns and Rows also swapped here.

	vtkSmartPointer<vtkCellArray> newOrderCells = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPoints> newOrderPoints = vtkSmartPointer<vtkPoints>::New();
	vtkIdType baseCell;
	vtkIdType basePoint;
	vtkIdType cellId;
	vtkIdType pointId;




	if (this->Rotations == 1) //smc
	{

		baseCell = this->Rows - 1;
		cellId =  baseCell;
		for (int i = 0; i < cells->GetNumberOfCells(); i++)
		{
			newOrderCells->InsertNextCell(subdividedQuad->GetCell(cellId));
			cellId += this->Rows;

			if ((i + 1) % (this->Columns) == 0)
			{
				baseCell--;
				cellId = baseCell;
			}
		}

		basePoint = this->Rows;
		pointId = basePoint;
		for (int i = 0; i < oldPoints->GetNumberOfPoints(); i++)
		{
			newOrderPoints->InsertNextPoint(oldPoints->GetPoint(pointId));
			pointId += (this->Rows + 1);

			if ((i + 1) % (this->Columns * 2 + 1) == 0)
			{
				basePoint--;
				pointId = basePoint;
			}

		}


	}
	else if (this->Rotations == 2)//ec TODO: doc, rows and columns reversed here as they are reversed earlier (only for ec)
	{
		baseCell = this->Columns * this->Rows;
		cellId = baseCell;
		for (int i = 0; i < cells->GetNumberOfCells(); i++)
		{
			newOrderCells->InsertNextCell(subdividedQuad->GetCell(cellId - 1));
			cellId--;
		}
	}

	else if (this->Rotations == 3)//smc
	{
		baseCell = this->Rows * (this->Columns - 1);
		cellId = this->Rows * (this->Columns - 1);
		for (int i = 0; i < cells->GetNumberOfCells(); i++)
		{
			newOrderCells->InsertNextCell(subdividedQuad->GetCell(cellId));

			if ((i + 1) % (this->Columns) == 0)
			{
				baseCell++;
				cellId = baseCell;
			}
			else
			{
				cellId -= this->Rows;
			}
		}
	}

	if (this->Rotations != 0)
	{
		//subdividedQuad->SetPolys(newOrderCells);
		//subdividedQuad->SetPoints(newOrderPoints);
	}

	output->ShallowCopy(subdividedQuad);
	return 1;
}

void vtkReorderSubdivideQuad::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Number of rows: " << this->Rows << "\n";
	os << indent << "Number of columns: " << this->Columns << "\n";
}
