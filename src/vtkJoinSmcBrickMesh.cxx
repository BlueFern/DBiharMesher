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
	this->Flat = true;
	this->AxialQuads = 0;
	this->CircQuads = 0;
	this->Branches = -1;
	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkJoinSmcBrickMesh::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
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

	for (int branchId = 0; branchId < this->Branches; branchId++)
	{
		quadId = 0;

		int branchOffset = branchId * cellsInQuad * quadsInBranch;
		fixes = this->Rows / 2; // Number of corrections required in a single quad.
		// First cell needed to be corrected in a quad, in each branch.
		int startPos = this->Columns + branchOffset;

		cell1  = startPos;
		cell2 = cell1 - cellsInQuad + this->Columns - 1;

		// Do each branch individually, the last iteration is most of the saddle
		while (quadId < quadsInBranch)
		{
			if (quadId % this->CircQuads == 0 && !this->Flat)
			{
				cell2 = cell1 + (cellsInQuad * (this->CircQuads - 1)) + this->Columns - 1;
			}

			vtkIdType newPoint[6];
			input->GetCellPoints(cell1, cell1Points);
			input->GetCellPoints(cell2, cell2Points);

			newPoint[0] = cell2Points->GetId(0);
			newPoint[1] = cell2Points->GetId(1);

			newPoint[2] = cell1Points->GetId(0);
			newPoint[3] = cell1Points->GetId(5);

			newPoint[4] = cell2Points->GetId(3);
			newPoint[5] = cell2Points->GetId(5);

			input->ReplaceCell(cell2, 6, newPoint);
			fixes--;
			// Finished in this quad, move to next.
			if (fixes == 0)
			{
				quadId++;
				cell1 = startPos + (cellsInQuad * quadId);
				cell2 = cell1 - cellsInQuad + this->Columns - 1;
				fixes = this->Rows / 2;
			}
			else
			{
				cell1 += this->Columns * 2;
				cell2 += this->Columns * 2;
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
