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
 * into a number of segments of that size (and all odd). As a result of this process
 * the sizes tend to be less than or equal to the partition length in all cases other
 * than spines over bifurcations and end points.
 *
 * Both straight sections and their reverses are added as cells. Bifurcations
 * are traversed by going down one branch (towards the bifurcation), and up the next.
 * This means there are as many segments/spines created as there are cells attached
 * to the bifurcation, and there will be duplication of points between cells.
 *
 * vtkPolyData is returned which has cells/lines as these segments, and points and
 * point data identical to the input.
 *
 */

#include <sstream>
#include <algorithm>
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
#include <vtkPointData.h>

#include "vtkCentrelinePartitioner.h"

vtkStandardNewMacro(vtkCentrelinePartitioner);
const int vtkCentrelinePartitioner::minEdgePoints = 5;

vtkCentrelinePartitioner::vtkCentrelinePartitioner()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
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
 * sized segments of the tree.
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

	vtkIdType branchId = 0;
	vtkIdType pointId = 0;
	vtkIdType localId = 0;

	cellIdList = input->GetCell(branchId)->GetPointIds();
	vtkIdType localEndPoint = cellIdList->GetNumberOfIds() - 1;
	vtkIdType cellSize = cellIdList->GetNumberOfIds();
	input->GetPointCells(cellIdList->GetId(localEndPoint), connectedCellIds);

	// If input vtkPolyData has just one line.
	if (connectedCellIds->GetNumberOfIds() == 1)
	{
		while (localId < cellSize)
		{
			pointId = cellIdList->GetId(localId);
			currentSegment->InsertUniqueId(pointId);
			localId++;
		}

		reverseIdList(currentSegment, reversedSpine);
		segments->InsertNextCell(currentSegment);
		segments->InsertNextCell(reversedSpine);
		output->SetLines(segments);

		return 1;
	}

	bool bifurcation = true;
	bool straightSegments = false;
	int actualLength = 0;
	int padding = 0;
	int branchStartingPoint[input->GetNumberOfCells()][2];
	int offset = 0;
	int numberOfSections = 0;

	int sections = cellSize / this->PartitionLength + (cellSize % this->PartitionLength !=  0); // Ceiling division.

	// Must have an odd number of sections for each section to be odd.
	if (sections % 2 == 0)
	{
		sections++;
	}

	int sectionLength = cellSize / sections;
	if (sectionLength % 2 == 0)
	{
		sectionLength++;
	}

	if (sections > 1)
	{
		straightSegments = true;
	}
	branchStartingPoint[0][0] = sectionLength;
	branchStartingPoint[0][1] = sections;

	input->GetPointCells(cellIdList->GetId(localEndPoint), connectedCellIds);
	int numberOfBranches = connectedCellIds->GetNumberOfIds();


	// Begin creating segments.

	for (int i = 1; i <= input->GetNumberOfCells(); i++)
	{
		currentSegment->Reset();
		joinedIdLists->Reset();
		reversedSpine->Reset();
		endSegment->Reset();

		if (straightSegments)
		{
			// If we're not on the trunk/inlet we need to begin from the end of the spine that
			// went down this branch already.

			actualLength = cellSize - 2 * branchStartingPoint[i-1][0]; // Reserved space at both ends of branch.
			if (!bifurcation || i == 1) // Only one reserved area to worry about.
			{
				actualLength = cellSize - branchStartingPoint[i-1][0];
			}

			sections = actualLength / this->PartitionLength + (actualLength % this->PartitionLength !=  0); // Ceiling division.

			// Must have an even number of sections for each section to be odd (duplicate points between segments).
			if (sections % 2 != 0)
			{
				sections++;
			}

			sectionLength = actualLength / sections;

			//exactDivide = (actualLength % sections == 0) && (sectionLength % 2 != 0);
			if (sectionLength % 2 == 0)
			{
				sectionLength++;
			}
			offset = 0;
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

				// Check if this branch ends in a bifurcation.
				input->GetPointCells(cellIdList->GetId(cellSize - 1), nextEndPoint);
				if (nextEndPoint->GetNumberOfIds() == 1)
				{
					bifurcation = false;
				}
				if (nextEndPoint->GetNumberOfIds() > 2)
				{
					bifurcation = true;
				}

				sections = cellSize / this->PartitionLength + (cellSize % this->PartitionLength !=  0);
				// With cells resampled with an odd number of points, must have at least 2 sections between two bifurcations.
				if (sections < 2)
				{
					sections = 2;
				}
				// Must have an even number of sections for each section to be odd (duplicate points between segments).
				else if (sections % 2 != 0)
				{
					sections++;
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
		if (i == input->GetNumberOfCells())
		{
			break;
		}

		cellIdList = input->GetCell(i)->GetPointIds();
		localEndPoint = cellIdList->GetNumberOfIds() - 1;

		// Does the next cell end in a bifurcation?
		input->GetPointCells(cellIdList->GetId(localEndPoint), connectedCellIds);
		numberOfBranches = connectedCellIds->GetNumberOfIds();

		if (numberOfBranches == 1)
		{
			bifurcation = false;
		}
		if (numberOfBranches > 2)
		{
			bifurcation = true;
		}

		cellSize = cellIdList->GetNumberOfIds();

		localId = branchStartingPoint[i][0] - 1;
		numberOfSections = branchStartingPoint[i][1];
		if (numberOfSections > 2 || (numberOfSections > 1 && bifurcation == false))
		{
			straightSegments = true;
		}
		else
		{
			straightSegments = false;
		}
		padding = branchStartingPoint[i][0] - 1;
	}

	output->SetLines(segments);

	return 1;
}

void vtkCentrelinePartitioner::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);

	//os << indent << "Minimum points for an edge: " << this->minEdgePoints << "\n";
	os << indent << "Given bound: " << this->PartitionLength << "\n";
}
