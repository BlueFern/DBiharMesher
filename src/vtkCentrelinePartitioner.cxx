#include <sstream>
#include <algorithm>
#include <vtkSmartPointer.h>
#include <vtkCallbackCommand.h>
#include <vtkIdTypeArray.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkPolyVertex.h>
#include <vtkVertex.h>
#include <vtkCellArray.h>
#include <vtkPriorityQueue.h>
#include <vtkPointData.h>

#include "vtkDbiharStatic.h"
#include "vtkCentrelinePartitioner.h"

vtkStandardNewMacro(vtkCentrelinePartitioner);
vtkCxxSetObjectMacro(vtkCentrelinePartitioner, EndPoints, vtkIdList);

vtkCentrelinePartitioner::vtkCentrelinePartitioner()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
	this->EndPoints = vtkIdList::New();
	this->PartitionLength = 0;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);

}

/**
 * Fills a second list (reversedSpine) with the IDs in the first (spine), in reverse order.
 */
void vtkCentrelinePartitioner::reverseIdList(vtkSmartPointer<vtkIdList> spine, vtkSmartPointer<vtkIdList> reversedSpine)
{
	reversedSpine->Reset(); // Ensure this list is initially empty.
	int spineId = spine->GetNumberOfIds() - 1;
	while (spineId >= 0)
	{
		reversedSpine->InsertNextId(spine->GetId(spineId));
		spineId--;
	}
}

/**
 * Fills the third empty list with the the first two joined together (second appended to the first).
 */
void vtkCentrelinePartitioner::joinIdLists(vtkSmartPointer<vtkIdList> previous, vtkSmartPointer<vtkIdList> current,
				   vtkSmartPointer<vtkIdList> joined)
{
	joined->Reset(); // Ensure the joined list is empty to begin.
	joined->DeepCopy(previous);
	for (int i = 0; i < current->GetNumberOfIds(); i++)
	{
		joined->InsertUniqueId(current->GetId(i));
	}
}

int vtkCentrelinePartitioner::RequestData(vtkInformation *vtkNotUsed(request),
										vtkInformationVector **inputVector,
										vtkInformationVector *outputVector)
{

	vtkPolyData *input = vtkPolyData::GetData(inputVector[0],0);
	input->BuildLinks();

	int numOfCells = input->GetNumberOfCells();
	vtkPolyData *output = vtkPolyData::GetData(outputVector,0);

	output->SetPoints(input->GetPoints());
	output->GetPointData()->SetScalars(input->GetPointData()->GetScalars());
	output->GetPointData()->SetVectors(input->GetPointData()->GetVectors());

	vtkSmartPointer<vtkCellArray> segments = vtkSmartPointer<vtkCellArray>::New();



	vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> endSegment = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> previousSegment = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> currentSegment = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> joinedIdLists = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> reversedSpine = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> connectedCellIds = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> tmpConnectedCellIds = vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkIdList> nextEndPoint = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> endPointsTmp = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> endPointsCopy = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> startingCell = vtkSmartPointer<vtkIdList>::New();

	// For setting vertices
	vtkSmartPointer<vtkCellArray> verticesArray = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPolyVertex> bifurcations = vtkSmartPointer<vtkPolyVertex>::New();
	vtkSmartPointer<vtkPolyVertex> endPoints = vtkSmartPointer<vtkPolyVertex>::New();

	if (this->EndPoints->GetNumberOfIds() == 0)
	{
		this->EndPoints->InsertNextId(0);
	}
	input->GetPointCells(this->EndPoints->GetId(0), startingCell);
	endPointsCopy->DeepCopy(this->EndPoints);

	if (startingCell->GetNumberOfIds() > 1)
	{
		vtkErrorMacro("Must not start on a bifurcation");
		exit(EXIT_FAILURE);
	}
	vtkIdType branchId = startingCell->GetId(0); // Start point. 0 if unspecified by user.

	cellIdList = input->GetCell(branchId)->GetPointIds();
	vtkIdType localId = vtkDbiharStatic::GetPosition(cellIdList, endPointsCopy->GetId(0));

	assert(localId != -1);
	endPointsCopy->DeleteId(endPointsCopy->GetId(0));

	vtkIdType localEndPoint = 0;

	endPointsTmp->DeepCopy(endPointsCopy);
	endPointsTmp->IntersectWith(cellIdList);
	if (endPointsTmp->GetNumberOfIds() == 1) // Still another end point on this cell.
	{
		localEndPoint = vtkDbiharStatic::GetPosition(cellIdList, endPointsTmp->GetId(0));
		endPointsCopy->DeleteId(endPointsTmp->GetId(0));
	}
	else
	{
		localEndPoint = cellIdList->GetNumberOfIds() - 1;
	}

	vtkIdType pointId = cellIdList->GetId(localId);
	endPoints->GetPointIds()->InsertNextId(pointId);

	vtkIdType cellSize = localEndPoint + 1 - localId;

	input->GetPointCells(cellIdList->GetId(localEndPoint), connectedCellIds);

	bool bifurcation = true;
	bool straightSegments = false;

	// No branches to traverse if last point it an endpoint or another specified endpoint exists in this cell.
	bool oneCell = (connectedCellIds->GetNumberOfIds() == 1 || endPointsTmp->GetNumberOfIds() == 1);

	int branchStartingPoint[input->GetNumberOfCells()][2];
	int offset = 0;
	int padding = localId;
	int numberOfSections = 0;
	int actualLengthCopy = 0;
	int sectionLength = 0;
	int stage = 1; // For rough progress tracking
	int sections = 0;
	int actualLength = cellSize - this->PartitionLength / 2;
	if (actualLength <= 0)
	{
		sections = 1;
	}
	else
	{
		sections = actualLength / this->PartitionLength + (actualLength % this->PartitionLength != 0);
	}


	// If only dealing with one cell we do not need to account for the space before a bifurcation
	// (this->PartitionLength / 2).
	if (oneCell)
	{
		sections = cellSize / this->PartitionLength + (cellSize % this->PartitionLength !=  0);
		straightSegments = true;
		bifurcation = false;

		// Endpoints to add to polyVertex.
		if (connectedCellIds->GetNumberOfIds() > 1)
		{
			endPoints->GetPointIds()->InsertNextId(endPointsTmp->GetId(0));
		}
		else
		{
			endPoints->GetPointIds()->InsertNextId(cellIdList->GetId(localEndPoint));
		}
	}

	if (sections == 1)
	{
		sectionLength = cellSize;
	}
	else
	{
		sectionLength = actualLength / sections;
	}

	if (sectionLength % 2 == 0)
	{
		sectionLength++;
	}

	straightSegments = sections > 1 || oneCell;

	branchStartingPoint[startingCell->GetId(0)][0] = this->PartitionLength / 2;
	branchStartingPoint[startingCell->GetId(0)][1] = sections;

	input->GetPointCells(cellIdList->GetId(localEndPoint), connectedCellIds);
	int numberOfBranches = connectedCellIds->GetNumberOfIds();
	vtkSmartPointer<vtkPriorityQueue> branchesToExplore = vtkSmartPointer<vtkPriorityQueue>::New();

	int priority = input->GetNumberOfCells();
	for (int i = 1; i < numberOfBranches; i++)
	{
		branchesToExplore->Insert(priority,(connectedCellIds->GetId(i)));
	}
	priority--;

	// Begin creating segments.

	while (true)
	{
		currentSegment->Reset();
		joinedIdLists->Reset();
		reversedSpine->Reset();
		endSegment->Reset();

		if (straightSegments)
		{
			// If we're not on the trunk/inlet we need to begin from the end of the spine that
			// went down this branch already.

			offset = 0;
			actualLength = cellSize - 2 * branchStartingPoint[branchId][0]; // Reserved space at both ends of branch.
			if (!bifurcation || branchId == startingCell->GetId(0)) // Only one reserved area to worry about.
			{
				if (oneCell)
				{
					actualLength = cellSize;
				}
				else
				{
					actualLength = cellSize - branchStartingPoint[branchId][0];
				}
			}

			// Number of sections NOT including reserved areas around bifurcations.
			sections = actualLength / this->PartitionLength + (actualLength % this->PartitionLength != 0);
			sectionLength = actualLength / sections;

			if (sectionLength % 2 == 0)
			{
				sectionLength++;
			}

			for (int section = 0; section < sections; section++)
			{
				// localId is decremented for each section so that they join, offset is used to correct for this.
				while (localId < padding + (sectionLength * (section + 1)) - offset)
				{
					pointId = cellIdList->GetId(localId);
					currentSegment->InsertUniqueId(pointId);
					localId++;
				}
				offset++;

				// Coming up to an end point (not a bifurcation), collect remaining points.
				if (section + 1 == sections && !bifurcation)
				{
					while (localId <= localEndPoint)
					{
						pointId = cellIdList->GetId(localId);
						currentSegment->InsertUniqueId(pointId);
						localId++;
					}
					pointId = cellIdList->GetId(localId);

				}

				localId--;

				reverseIdList(currentSegment, reversedSpine);
				segments->InsertNextCell(currentSegment);
				segments->InsertNextCell(reversedSpine);
				currentSegment->Reset();
				reversedSpine->Reset();
			}
		}

		if (bifurcation)
		{
			// Collect points up to bifurcation.
			while (localId < localEndPoint)
			{
				pointId = cellIdList->GetId(localId);
				endSegment->InsertUniqueId(pointId);
				localId++;
			}
			pointId = cellIdList->GetId(localId);
			bifurcations->GetPointIds()->InsertNextId(pointId);

			previousSegment->DeepCopy(endSegment);

			// Create spines around the bifurcation.
			// Assumes the branches are in order (for more than 2 branches).
			for (int j = 1; j < numberOfBranches; j++)
			{
				localId = 0;
				branchId = connectedCellIds->GetId(j);
				cellIdList = input->GetCell(branchId)->GetPointIds();
				cellSize = cellIdList->GetNumberOfIds();

				endPointsTmp->Reset();
				endPointsTmp->DeepCopy(endPointsCopy);
				endPointsTmp->IntersectWith(cellIdList);
				if (endPointsTmp->GetNumberOfIds() > 0)
				{
					if (endPointsTmp->GetNumberOfIds() > 1)
					{
						vtkWarningMacro("Multiple end points on same branch");
					}

					cellSize = endPointsTmp->GetId(0) - cellIdList->GetId(1) + 2;
				}

				sections = (cellSize + this->PartitionLength / 2) / this->PartitionLength + (cellSize % this->PartitionLength !=  0); // Ceiling division.

				// Does this branch end in a natural endpoint (no bifurcation)?
				localEndPoint = cellIdList->GetNumberOfIds() - 1;
				input->GetPointCells(cellIdList->GetId(localEndPoint), tmpConnectedCellIds);
				int endBranches = tmpConnectedCellIds->GetNumberOfIds();


				// Must have at least two sections between bifurcations (unless there is an endpoint).
				// Endpoint could either be a specified endpoint or natural one.
				if (sections == 1)
				{
					// First check for endpoints to add to polyVertex.

					if (endPointsTmp->GetNumberOfIds() > 0)
					{
						endPoints->GetPointIds()->InsertNextId(endPointsTmp->GetId(0));
					}
					else if (endBranches == 1)
					{
						endPoints->GetPointIds()->InsertNextId(cellIdList->GetId(localEndPoint));
					}

					else
					{
						sections = 2;
					}
					sectionLength = cellSize / sections;
				}

				else
				{
					sectionLength = this->PartitionLength / 2;
				}

				if (sectionLength % 2 == 0)
				{
					sectionLength++;
				}

				branchStartingPoint[branchId][0] = sectionLength;
				branchStartingPoint[branchId][1] = sections;

				while (localId < sectionLength)
				{
					pointId = cellIdList->GetId(localId);
					currentSegment->InsertUniqueId(pointId);
					localId++;

					// If the end point is in the spine around the bifurcation, we do not need to explore it later.
					if (localId == cellSize)
					{
						branchesToExplore->DeleteId(branchId);
						endPointsCopy->DeleteId(pointId);
					}
				}

				// In the first iteration the points are already in ascending order.
				if (j == 1)
				{
					joinIdLists(previousSegment, currentSegment, joinedIdLists);
				}
				else
				{
					reverseIdList(previousSegment,reversedSpine);
					joinIdLists(reversedSpine, currentSegment, joinedIdLists);
				}
				segments->InsertNextCell(joinedIdLists);
				previousSegment->Reset();
				previousSegment->DeepCopy(currentSegment);
				currentSegment->Reset();
				joinedIdLists->Reset();
				reversedSpine->Reset();
			}

			// Final spine around bifurcation.
			joinIdLists(endSegment, previousSegment, joinedIdLists); // Is moving the opposite direction.
			reverseIdList(joinedIdLists, reversedSpine);
			segments->InsertNextCell(reversedSpine);
			currentSegment->Reset();
			joinedIdLists->Reset();
			previousSegment->Reset();
			reversedSpine->Reset();
			endSegment->Reset();
		}

		// Get ready for next iteration.

		// No data to get during last iteration. Break early.
		if (branchesToExplore->GetNumberOfItems() == 0)
		{
			if (endPointsCopy->GetNumberOfIds() != 0)
			{
				vtkWarningMacro("One or more end points were not reached.");
			}
			break;
		}

		branchId = branchesToExplore->Pop();
		cellIdList = input->GetCell(branchId)->GetPointIds();
		endPointsTmp->Reset();
		endPointsTmp->DeepCopy(endPointsCopy);
		endPointsTmp->IntersectWith(cellIdList);


		if (endPointsTmp->GetNumberOfIds() > 0)
		{
			// Next branch has an endpoint down it - but past the reserved area around the bifurcation
			localEndPoint = endPointsTmp->GetId(0) - cellIdList->GetId(1) + 1;
			endPoints->GetPointIds()->InsertNextId(cellIdList->GetId(localEndPoint));

			endPointsCopy->DeleteId(endPointsTmp->GetId(0));

			bifurcation = false;
		}
		else
		{
			localEndPoint = cellIdList->GetNumberOfIds() - 1;


			// Does the next cell end in a bifurcation?
			input->GetPointCells(cellIdList->GetId(localEndPoint), connectedCellIds);
			numberOfBranches = connectedCellIds->GetNumberOfIds();

			if (numberOfBranches == 1)
			{
				endPoints->GetPointIds()->InsertNextId(cellIdList->GetId(localEndPoint));

			}

			for (int i = 1; i < numberOfBranches; i++) // i = 0 is the previous branch.
			{
				branchesToExplore->Insert(priority,(connectedCellIds->GetId(i)));
			}
			priority--;

			bifurcation = numberOfBranches > 2;

		}

		cellSize = localEndPoint + 1;

		localId = branchStartingPoint[branchId][0] - 1;
		numberOfSections = branchStartingPoint[branchId][1];

		straightSegments = (numberOfSections > 2 || (numberOfSections > 1 && bifurcation == false));

		padding = branchStartingPoint[branchId][0] - 1;

		this->UpdateProgress(static_cast<double>(stage++) / static_cast<double>(numOfCells));
	}


	verticesArray->InsertNextCell(endPoints);
	if (bifurcations->GetNumberOfPoints() > 0)
	{
		verticesArray->InsertNextCell(bifurcations);
	}

	output->SetVerts(verticesArray);
	output->SetLines(segments);

	return 1;
}

void vtkCentrelinePartitioner::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkCentrelinePartitioner* filter = static_cast<vtkCentrelinePartitioner *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}

void vtkCentrelinePartitioner::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Given bound: " << this->PartitionLength << "\n";
	int numIds = this->EndPoints->GetNumberOfIds();
	int endPoints[numIds];
	os << indent << "End Points: ";
	for (int i = 0; i < numIds; i++)
	{
		os << (this->EndPoints->GetId(i));
		if (i + 1 < numIds)
		{
			os << ", ";
		}
	}
	os << ".";
}
