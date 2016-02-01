#include <map>

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkInformationVector.h>
#include <vtkIdList.h>
#include <vtkCallbackCommand.h>

#include "vtkJoinSmcBrickMesh.h"

vtkStandardNewMacro(vtkJoinSmcBrickMesh);

vtkJoinSmcBrickMesh::vtkJoinSmcBrickMesh()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Columns = 0;
	this->Rows = 0;
	this->Flat = false;
	this->AxialQuads = 0;
	this->CircQuads = 0;
	this->Branches = 0;
	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkJoinSmcBrickMesh::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	// Basic testing of input parameters.
	if (this->Columns <= 0 || this->Rows <= 0)
	{
		vtkErrorMacro("Both Columns and Rows must be positive numbers.");
		exit(EXIT_FAILURE);
	}
	if (this->AxialQuads <= 0 || this->CircQuads <= 0)
	{
		vtkErrorMacro("Both the number of circumferential and axial quads must be positive numbers.");
		exit(EXIT_FAILURE);
	}
	if (this->Branches == 0)
	{
		vtkErrorMacro("Must specify the number of branches.");
		exit(EXIT_FAILURE);
	}

	vtkSmartPointer<vtkIdList> cell1Points = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> cell2Points = vtkSmartPointer<vtkIdList>::New();
	int cellsInQuad = this->Rows * this->Columns;
	int quadsInBranch = this->CircQuads * this->AxialQuads;
	vtkIdType cell1;
	vtkIdType cell2;
	vtkIdType quadId;
	int fixes;
	int startPos;
	int branchOffset;

	for (int branchId = 0; branchId < this->Branches; branchId++)
	{
		quadId = 0;

		branchOffset = branchId * cellsInQuad * quadsInBranch;
		fixes = this->Rows / 2; // Number of corrections/joins required in a single quad.

		startPos = branchOffset;

		cell1  = startPos + cellsInQuad;
		cell2 = cell1 - cellsInQuad + this->Columns - 1;

		while (quadId < quadsInBranch)
		{
			vtkIdType newPoint[6];
			input->GetCellPoints(cell1, cell1Points);
			input->GetCellPoints(cell2, cell2Points);

			// The quads in the top half for each branch have a different pattern to those in the bottom half.
			// The connectivity if therefore different, and they must be done seperately.

			if (quadId % this->CircQuads < this->CircQuads / 2) // Top half.
			{
				newPoint[0] = cell2Points->GetId(0);
				newPoint[1] = cell2Points->GetId(1);

				newPoint[2] = cell1Points->GetId(0);
				newPoint[3] = cell1Points->GetId(3);
				newPoint[4] = cell2Points->GetId(3);
				newPoint[5] = cell2Points->GetId(5);

				input->ReplaceCell(cell2, 6, newPoint);

			}
			else if (quadId % this->CircQuads >= this->CircQuads / 2) // Bottom half.
			{
				newPoint[0] = cell2Points->GetId(1);
				newPoint[1] = cell1Points->GetId(0);
				newPoint[2] = cell1Points->GetId(2);
				newPoint[3] = cell1Points->GetId(3);
				newPoint[4] = cell1Points->GetId(4);
				newPoint[5] = cell2Points->GetId(2);

				input->ReplaceCell(cell1, 6, newPoint);
			}

			fixes--;
			// Finished in this quad, move to next.
			if (fixes == 0)
			{
				quadId++;
				fixes = this->Rows / 2;

				// Nothing exists to be fixed on the last row (given fixing takes place for a quad with its neighbour
				// above it). We move to the next quad.
				if  ((quadId + 1) % this->CircQuads == 0 )
				{
					quadId++;
					cell1 += cellsInQuad + this->Columns;
					cell2 += cellsInQuad + this->Columns;
				}
				else
				{
					cell1 = startPos + (cellsInQuad * (quadId + 1));
					cell2 = cell1 - cellsInQuad + this->Columns - 1;
				}


				// The bottom half quads (circumferentially) have a different pattern - the cells
				// to be joined are in different positions.
				if ((quadId + 1) % this->CircQuads >= this->CircQuads / 2)
				{
					cell1 += this->Columns;
					cell2 += this->Columns;

				}

				// The existance of 'short' cells along this boundary makes no joining of cells necessary.
				// We simply skip the middle row.
				if ((quadId + 1) % this->CircQuads == this->CircQuads / 2)
				{
					quadId++;
					cell1 += cellsInQuad;
					cell2 += cellsInQuad;
				}
			}
			else
			{
				cell1 += this->Columns * 2;
				cell2 += this->Columns * 2;
			}
		}
	}

	// If the mesh is not flat, the two different patterns in the top and bottom halves of each branch
	// result in gaps. To correct this cells along the join are stretched.
	if (!this->Flat)
	{
		vtkIdType newPoint[4];
		quadId = 0;

		for (int branchId = 0; branchId < this->Branches; branchId++)
		{
			int ringId = 0;
			branchOffset = branchId * cellsInQuad * quadsInBranch;
			startPos = branchOffset;
			cell1 = startPos;
			cell2 = cell1 + cellsInQuad * (this->CircQuads - 1) + this->Columns - 1;
			fixes = this->Rows;


			while (ringId < this->AxialQuads)
			{
				vtkIdType newPoint[6];
				input->GetCellPoints(cell1, cell1Points);
				input->GetCellPoints(cell2, cell2Points);

				// Along this join, every second cell from the upper quad needs to be streched down.
				// Similarly, every other cell from the lower quad needs to be streched upwards.

				if (quadId % 2 == 0) // Stretched up.
				{

					newPoint[0] = cell2Points->GetId(1);
					newPoint[1] = cell1Points->GetId(1);
					newPoint[2] = cell1Points->GetId(2);
					newPoint[3] = cell2Points->GetId(2);

					input->ReplaceCell(cell1, 4, newPoint);
				}
				else // Stretched down.
				{
					newPoint[0] = cell2Points->GetId(0);
					newPoint[1] = cell1Points->GetId(0);
					newPoint[2] = cell1Points->GetId(3);
					newPoint[3] = cell2Points->GetId(3);

					input->ReplaceCell(cell2, 4, newPoint);
				}

				quadId++;
				fixes--;
				// Finished in this quad, move to next.
				if (fixes == 0)
				{
					ringId++;
					quadId = 0;
					cell1 = startPos + (ringId * cellsInQuad * this->CircQuads);
					cell2 = cell1 + cellsInQuad * (this->CircQuads - 1) + this->Columns - 1;
					fixes = this->Rows;
				}
				else
				{
					cell1 += this->Columns;
					cell2 += this->Columns;
				}
			}
		}
	}

	output->ShallowCopy(input);
	return 1;
}

void vtkJoinSmcBrickMesh::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Columns: " << this->Columns << "\n";
	os << indent << "Rows: " << this->Rows << "\n";
}

void vtkJoinSmcBrickMesh::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkJoinSmcBrickMesh* filter = static_cast<vtkJoinSmcBrickMesh *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}
