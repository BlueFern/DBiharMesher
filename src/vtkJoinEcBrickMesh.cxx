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
	this->Branches = 0;
	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkJoinEcBrickMesh::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
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

	output->DeepCopy(input);

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

		startPos = this->Columns * (this->Rows - 1) + branchOffset + 1;

		// The second and third branches (for such a mesh) have a different cell pattern to allow for a
		// correct tessellation. The start position is then different.
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

			// We do not correct cells around the sadle or at the very ends of daughter branches -
			// they are special cases handled later.
			if (quadId >= quadsInBranch - this->CircQuads)
			{
				break;
			}

			vtkIdType newPoint[6];
			output->GetCellPoints(cell1, cell1Points);
			output->GetCellPoints(cell2, cell2Points);


			if (branchId > 0) // Daughter branches.
			{
				newPoint[0] = cell1Points->GetId(3);
				newPoint[1] = cell1Points->GetId(2);

				newPoint[2] = cell2Points->GetId(2);
				newPoint[3] = cell2Points->GetId(4);

				newPoint[4] = cell2Points->GetId(5);
				newPoint[5] = cell2Points->GetId(0);
				output->ReplaceCell(cell2, 6, newPoint);
			}
			else // Trunk.
			{
				newPoint[0] = cell1Points->GetId(0);
				newPoint[1] = cell1Points->GetId(1);

				newPoint[2] = cell1Points->GetId(2);
				newPoint[3] = cell2Points->GetId(1);

				newPoint[4] = cell2Points->GetId(0);
				newPoint[5] = cell1Points->GetId(5);
				output->ReplaceCell(cell1, 6, newPoint);
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
