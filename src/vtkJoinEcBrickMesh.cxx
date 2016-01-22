#include <map>

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkInformationVector.h>
#include <vtkIdList.h>
#include <vtkCallbackCommand.h>

#include "vtkJoinEcBrickMesh.h"

vtkStandardNewMacro(vtkJoinEcBrickMesh);

vtkJoinEcBrickMesh::vtkJoinEcBrickMesh()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Columns = 0;
	this->Rows = 0;
	this->AxialQuads = 0;
	this->CircQuads = 0;
	this->Branches = 3;
	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkJoinEcBrickMesh::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	vtkSmartPointer<vtkIdList> cell1Points = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> cell2Points = vtkSmartPointer<vtkIdList>::New();

	int cellsInQuad = this->Rows * this->Columns;
	int quadsInBranch = this->CircQuads * this->AxialQuads;
	vtkIdType cell1;
	vtkIdType cell2;
	vtkIdType quadId;
	int fixes;
	int startPos;

	for (int branchId = 0; branchId < this->Branches; branchId++)
	{
		int branchOffset = branchId * cellsInQuad * quadsInBranch;
		fixes = this->Columns / 2; // Number of corrections required in a single quad.
		// First cell needed to be corrected in a quad, in each branch.

		//doc TODO
		startPos = this->Columns * (this->Rows - 1) + branchOffset + 1;
		if (branchId > 0)
		{
			startPos--;
		}

		cell1  = startPos;
		cell2 = cell1 + this->CircQuads * cellsInQuad - (this->Columns * (this->Rows - 1));

		quadId = 0;

		// Do each branch individually
		while (quadId < quadsInBranch)
		{

			if (quadId >= quadsInBranch - this->CircQuads)
			{
				break;
			}
			vtkIdType newPoint[6];
			input->GetCellPoints(cell1, cell1Points);
			input->GetCellPoints(cell2, cell2Points);

			if (branchId > 0)
			{
				newPoint[0] = cell1Points->GetId(5);
				newPoint[1] = cell2Points->GetId(3);

				newPoint[2] = cell2Points->GetId(5);
				newPoint[3] = cell2Points->GetId(0);

				newPoint[4] = cell2Points->GetId(1);
				newPoint[5] = cell1Points->GetId(0);
				input->ReplaceCell(cell2, 6, newPoint);
			}
			else
			{
				newPoint[0] = cell1Points->GetId(0);
				newPoint[1] = cell1Points->GetId(1);

				newPoint[2] = cell2Points->GetId(0);
				newPoint[3] = cell2Points->GetId(5);

				newPoint[4] = cell1Points->GetId(4);
				newPoint[5] = cell1Points->GetId(5);
				input->ReplaceCell(cell1, 6, newPoint);
			}

			fixes--;

			// Finished in this quad, move to next.
			if (fixes == 0)
			{
				quadId++;
				cell1 = startPos + (cellsInQuad * quadId);
				cell2 = (cellsInQuad * quadId) + (cellsInQuad * this->CircQuads)+ branchOffset + 1;
				if (branchId > 0)
				{
					cell2--;
				}
				fixes = this->Columns / 2;
			}
			else
			{
				cell1 += 2;
				cell2 += 2;
			}
		}
	}
#if 0
	// Fix the remaining part of the saddle by stretching cells across.
	// Work from both sides at once.

	vtkIdType baseQuad1 = this->CircQuads * this->AxialQuads + this->CircQuads / 2;
	vtkIdType baseQuad2 = this->CircQuads * this->AxialQuads * 2 + this->CircQuads / 2 - 1;

	vtkIdType baseCell1 = baseQuad1 * cellsInQuad;
	vtkIdType baseCell2 = baseQuad2 * cellsInQuad + this->Columns - 1;

	cell1 = baseCell1;
	cell2 = baseCell2;

	fixes = this->Columns;
	quadId = 0;

	while (quadId < this->CircQuads / 2)
	{
		vtkIdType newCell[6];

		input->GetCellPoints(cell1, cell1Points);
		input->GetCellPoints(cell2, cell2Points);

		if (fixes % 2 == 0)
		{
			newCell[0] = cell2Points->GetId(5);
			newCell[1] = cell1Points->GetId(0);
			newCell[2] = cell1Points->GetId(1);
			newCell[3] = cell1Points->GetId(3);
			newCell[4] = cell1Points->GetId(5);
			newCell[5] = cell2Points->GetId(0);

			input->ReplaceCell(cell1, 6, newCell);
		}
		else
		{
			newCell[0] = cell2Points->GetId(3);
			newCell[1] = cell2Points->GetId(5);
			newCell[2] = cell1Points->GetId(0);
			newCell[3] = cell1Points->GetId(5);
			newCell[4] = cell2Points->GetId(0);
			newCell[5] = cell2Points->GetId(1);

			input->ReplaceCell(cell2, 6, newCell);
		}

		fixes--;
		// Finished in this quad, move to next.
		if (fixes == 0)
		{
			quadId++;
			cell1 = baseCell1 + (cellsInQuad * quadId);
			cell2 = baseCell2 - (cellsInQuad * quadId);
			fixes = this->Columns;
		}
		else
		{
			cell1 += 1;
			cell2 -= 1;
		}

	}
#endif
	output->ShallowCopy(input);
	return 1;
}

void vtkJoinEcBrickMesh::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Columns: " << this->Columns << "\n";
	os << indent << "Rows: " << this->Rows << "\n";
}

void vtkJoinEcBrickMesh::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkJoinEcBrickMesh* filter = static_cast<vtkJoinEcBrickMesh *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}
