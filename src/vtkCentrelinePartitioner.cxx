/**
 * This filter partitions given vtkPolyData into segments using a user set partition
 * length as a guide. This bound should represent how many points (roughly) are to be
 * included to or from a bifurcation point to avoid having non-smooth boundaries between
 * branches. Therefore the segment size for a spine over such a bifurcation will be
 * closer to twice the partition length.
 *
 * Less important are the straight segments at the ends or between bifurcations. Their
 * size will be roughly the input bound, and always greater than the minimum number
 * of points the dbihar patch filter requires to build an edge (so long as the input
 * data has enough points in each cell).
 *
 * The size of the segments is only ever rough due to the requirement for each segment
 * to have odd length. Using the given bound, the program attempts to divide a cell/branch
 * into a number of segments of that size (odd length). As a result of this process
 * the sizes tend to be less than or equal to the partition length in all cases other
 * than spines over bifurcations and end points.
 *
 * Both straight sections and their reverses are added as cells. Bifurcations
 * are traversed by going down one branch (towards the bifurcation), and up the next.
 * This means there are as many segments/spines created as there are cells attached
 * to the bifurcation, and there will be duplication of points between cells.
 *
 * EndPoints is an optional Id list that specifies points to build segments between.
 * The first point is where the partitioner should start, and all others are used as new
 * end points for the associated cell. The user is responsible for giving sensible end
 * points so that the cell sizes are still all odd. If points are given that belong
 * in cells that branched out before the given starting point they will be ignored.
 *
 * vtkPolyData is returned which has cells/lines as these segments, and points and
 * point data identical to the input.
 *
 */

#include <sstream>
#include <algorithm>
#include <string>
#include <vtkSmartPointer.h>
#include <vtkIdTypeArray.h>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkCellArray.h>
#include <vtkPriorityQueue.h>
#include <vtkPointData.h>

#include "vtkCentrelinePartitioner.h"

vtkStandardNewMacro(vtkCentrelinePartitioner);

vtkCentrelinePartitioner::vtkCentrelinePartitioner()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
}

/**
 * Finds the local Id of a global Id and the cell it belongs to.
 */
vtkIdType vtkCentrelinePartitioner::findLocalId(vtkSmartPointer<vtkIdList> list, vtkIdType pointId)
{
	for (vtkIdType id = 0; id < list->GetNumberOfIds(); id++)
	{
		if (list->GetId(id) == pointId)
		{
			return id;
		}
	}
	return -1;
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

/**
 * Main logic of the filter. Iteratively traverses the input centreline data and builds up cells that are roughly
 * similarly sized segments.
 */
int vtkCentrelinePartitioner::RequestData(vtkInformation *vtkNotUsed(request),
										vtkInformationVector **inputVector,
										vtkInformationVector *outputVector)
{
	vtkPolyData *input = vtkPolyData::GetData(inputVector[0],0);
	input->BuildLinks();

	vtkPolyData *output = vtkPolyData::GetData(outputVector,0);

	output->SetPoints(input->GetPoints());
	output->GetPointData()->SetScalars(input->GetPointData()->GetScalars());

	vtkSmartPointer<vtkCellArray> segments = vtkSmartPointer<vtkCellArray>::New();

	vtkSmartPointer<vtkIdList> cellIdList = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> endSegment = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> previousSegment = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> currentSegment = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> joinedIdLists = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> reversedSpine = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> connectedCellIds = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> nextEndPoint = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> endPointsTmp = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> endPointsCopy = vtkSmartPointer<vtkIdList>::New(); // For printing after endpoints are popped off.
	vtkSmartPointer<vtkIdList> startingCell = vtkSmartPointer<vtkIdList>::New();

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
	vtkIdType localId = findLocalId(cellIdList, endPointsCopy->GetId(0));
	assert(localId != -1);
	endPointsCopy->DeleteId(endPointsCopy->GetId(0));

	vtkIdType localEndPoint = 0;

	endPointsTmp->DeepCopy(endPointsCopy);
	endPointsTmp->IntersectWith(cellIdList);
	if (endPointsTmp->GetNumberOfIds() == 1) // Still another end point on this cell.
	{
		localEndPoint = findLocalId(cellIdList, endPointsTmp->GetId(0));
		endPointsCopy->DeleteId(endPointsTmp->GetId(0));
	}
	else
	{
		localEndPoint = cellIdList->GetNumberOfIds() - 1;
	}

	vtkIdType pointId = 0;
	vtkIdType cellSize = localEndPoint + 1 - localId;

	input->GetPointCells(cellIdList->GetId(localEndPoint), connectedCellIds);

	bool bifurcation = true;
	bool straightSegments = false;

	// Already deleted starting point, number of Ids is 1 for single cell partition.
	bool oneCell = (connectedCellIds->GetNumberOfIds() == 1 || endPointsTmp->GetNumberOfIds() == 1);
	int actualLength = 0;
	int padding = 0;
	int branchStartingPoint[input->GetNumberOfCells()][2];
	int offset = 0;
	int numberOfSections = 0;

	int sections = cellSize / this->PartitionLength + (cellSize % this->PartitionLength !=  0); // Ceiling division.

	int sectionLength = cellSize / sections;
	if (sectionLength % 2 == 0)
	{
		sectionLength++;
	}

	straightSegments = sections > 1;

	if (oneCell)
	{
		straightSegments = true;
		bifurcation = false;
	}

	branchStartingPoint[startingCell->GetId(0)][0] = sectionLength;
	branchStartingPoint[startingCell->GetId(0)][1] = sections;

	input->GetPointCells(cellIdList->GetId(localEndPoint), connectedCellIds);
	int numberOfBranches = connectedCellIds->GetNumberOfIds();
	vtkSmartPointer<vtkPriorityQueue> branchesToExplore = vtkSmartPointer<vtkPriorityQueue>::New();

	int priority = input->GetNumberOfCells();
	for (int i = 1; i < numberOfBranches; i++) // First item is backwards.
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
					offset -= localId;
				}
				else
				{
					actualLength = cellSize - branchStartingPoint[branchId][0];
				}
			}

			sections = actualLength / this->PartitionLength + (actualLength % this->PartitionLength !=  0); // Ceiling division.
			sectionLength = actualLength / sections;

			if (sectionLength % 2 == 0)
			{
				sectionLength++;
			}

			for (int section = 0; section < sections; section++)
			{
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

				sections = cellSize / this->PartitionLength + (cellSize % this->PartitionLength !=  0);

				// Must have at least two sections between bifurcations.
				if (sections < 2)
				{
					sections = 2;
				}

				sectionLength = cellSize / sections;
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
			localEndPoint = endPointsTmp->GetId(0) - cellIdList->GetId(1) + 1;
			endPointsCopy->DeleteId(endPointsTmp->GetId(0));
		}
		else
		{
			localEndPoint = cellIdList->GetNumberOfIds() - 1;
		}

		// Does the next cell end in a bifurcation?
		input->GetPointCells(cellIdList->GetId(localEndPoint), connectedCellIds);
		numberOfBranches = connectedCellIds->GetNumberOfIds();

		for (int i = 1; i < numberOfBranches; i++) // First item is backwards.
		{
			branchesToExplore->Insert(priority,(connectedCellIds->GetId(i)));
		}
		priority--;

		bifurcation = numberOfBranches > 2;

		cellSize = localEndPoint + 1;

		localId = branchStartingPoint[branchId][0] - 1;
		numberOfSections = branchStartingPoint[branchId][1];

		straightSegments = (numberOfSections > 2 || (numberOfSections > 1 && bifurcation == false));

		padding = branchStartingPoint[branchId][0] - 1;
	}

	output->SetLines(segments);

	return 1;
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
