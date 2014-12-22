/**
 * This filter partitions given vtkPolyData into segments with maximum size of
 * a given bound + a small fixed distance (so that straight segments are at least
 * 4 points long - a minimum requirement for the dbihar patch generator).
 *
 * Both straight sections and their reverses are added as cells. Bifurcations
 * are traversed by going down one branch (towards the bifurcation), and up the next.
 * This means there are as many segments/spines created as there are cells attached
 * to the bifurcation, and there will be duplication of points between cells.
 *
 * vtkPolyData is returned which has cells/lines as these segments, and points and
 * point data identical to the input.
 *
 *
 * Created by: Stewart Dowding
 * Version 1.1	 22/12/2014
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
#include <vtkCentrelinePartitioner.h>
#include <vtkPointData.h>

#include "showPolyData.h"

vtkStandardNewMacro(vtkCentrelinePartitioner);

vtkCentrelinePartitioner::vtkCentrelinePartitioner()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
}

int GetBound(bool bifurcation, int cellSize, int Bound)
{
	int actualBound = 0;
	int minEdgePoints = 4;
	int numberOfBifurcations = 1; //Need twice the space for two bifurcations.

	if (bifurcation)
	{
		numberOfBifurcations = 2;
	}

	if (cellSize <= numberOfBifurcations * Bound)
	{
		actualBound = cellSize / numberOfBifurcations;
	}
	else
	{
		actualBound = Bound;

		// Big enough for Bounds, but not leaving enough room for a straight segment.
		if (cellSize < numberOfBifurcations * Bound + minEdgePoints)
		{
			actualBound = cellSize / numberOfBifurcations;
		}
	}

	return actualBound;

}

void reverseIdList(vtkSmartPointer<vtkIdList> spine, vtkSmartPointer<vtkIdList> reversedSpine)
{
	int spineLength = spine->GetNumberOfIds() - 1;
	while (spineLength >= 0)
	{
		reversedSpine->InsertNextId(spine->GetId(spineLength));
		spineLength--;
	}
}

void joinIdLists(vtkSmartPointer<vtkIdList> previous, vtkSmartPointer<vtkIdList> current,
				   vtkSmartPointer<vtkIdList> joined)
{
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

	vtkIdType branchId;
	vtkIdType pointId = 0;
	vtkIdType localId = 0;

	cellIdList = input->GetCell(0)->GetPointIds();
	vtkIdType localEndPoint = cellIdList->GetNumberOfIds() - 1;
	vtkIdType endPoint = cellIdList->GetId(localEndPoint);
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

		segments->InsertNextCell(currentSegment);
		output->SetLines(segments);

		return 1;

	}

	bool bifurcation = true;
	bool straightSegments = false;
	int minEdgepoints = 4;
	int actualBound;
	int sections;
	int numberOfBranches = 0;

	actualBound = GetBound(false, cellSize, this->Bound);

	if (cellSize > 2 * this->Bound + minEdgepoints || cellSize == 2 * this->Bound)
	{
		straightSegments = true;
	}

	//-----------------------------------------------------------------------------------------------
	// Begin creating segments.
	//-----------------------------------------------------------------------------------------------

	for (int i = 1; i <= input->GetNumberOfCells(); i++)
	{
		if (straightSegments)
		{
			// How many full straight sections we can make before hitting 'bifurcation-reserved' area?
			sections = (cellSize - actualBound) / actualBound;

			for (int k = 0; k <= sections; k++)
			{
				while (localId < actualBound * (k + 1))
				{
					pointId = cellIdList->GetId(localId);
					currentSegment->InsertUniqueId(pointId);
					localId++;
				}

				localId--; // For connectivity in next cell.

				// Not pretty. Skips over a straight segment that was part of a previous bifurcation segment.
				if (currentSegment->GetNumberOfIds() > 0)
				{
					reverseIdList(currentSegment, reversedSpine);
					segments->InsertNextCell(currentSegment);
					segments->InsertNextCell(reversedSpine);
					currentSegment->Reset();
					reversedSpine->Reset();
				}
			}

			if (bifurcation)
			{
				// Add left over bit until hitting 'bifurcation-reserved' area.
				while (localId < cellSize - actualBound)
				{
					pointId = cellIdList->GetId(localId);
					currentSegment->InsertUniqueId(pointId);
					localId++;
				}

				if (currentSegment->GetNumberOfIds() > 0)
				{
					reverseIdList(currentSegment, reversedSpine);
					segments->InsertNextCell(currentSegment);
					segments->InsertNextCell(reversedSpine);
					currentSegment->Reset();
					reversedSpine->Reset();
				}
			}
		}

		//-----------------------------------------------------------------------------------------------
		// Collect points up to bifurcation/end point.
		//-----------------------------------------------------------------------------------------------

		while (localId <= localEndPoint)
		{
			pointId = cellIdList->GetId(localId);
			endSegment->InsertUniqueId(pointId);
			localId++;
		}

		input->GetPointCells(cellIdList->GetId(localEndPoint), connectedCellIds);
		numberOfBranches = connectedCellIds->GetNumberOfIds();

		if(numberOfBranches == 1)
		{
			bifurcation = false;
		}
		if (numberOfBranches > 2)
		{
			bifurcation = true;
		}

		if (bifurcation)
		{
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
				if(nextEndPoint->GetNumberOfIds() == 1)
				{
					bifurcation = false;
				}
				if (nextEndPoint->GetNumberOfIds() > 2)
				{
					bifurcation = true;
				}

				actualBound = GetBound(bifurcation, cellSize, this->Bound);

				while (localId < actualBound)
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

			joinIdLists(endSegment, previousSegment, joinedIdLists); // Is moving the opposite direction.
			reverseIdList(joinedIdLists, reversedSpine);
			segments->InsertNextCell(reversedSpine);
			currentSegment->Reset();
			joinedIdLists->Reset();
			previousSegment->Reset();
			reversedSpine->Reset();
			endSegment->Reset();
		}

		// Wasn't a bifurcation. Add the end segment and its reverse to the cell array.
		else
		{
			segments->InsertNextCell(endSegment);
			reverseIdList(endSegment, reversedSpine);
			segments->InsertNextCell(reversedSpine);
			endSegment->Reset();
			reversedSpine->Reset();
		}

		//------------------------------------------------------------------------------------------------
		// Get ready for next iteration.
		//-----------------------------------------------------------------------------------------------

		// No data to get during last iteration. Break early.
		// For-loop format necessary because of the difference with the trunk cell.
		if (i == input->GetNumberOfCells())
		{
			break;
		}

		cellIdList = input->GetCell(i)->GetPointIds();
		localEndPoint = cellIdList->GetNumberOfIds() - 1;
		endPoint = cellIdList->GetId(localEndPoint);
		input->GetPointCells(endPoint, connectedCellIds);

		// Does the next cell end in a bifurcation?
		if(connectedCellIds->GetNumberOfIds() == 1)
		{
			bifurcation = false;
		}
		if (connectedCellIds->GetNumberOfIds() > 2)
		{
			bifurcation = true;
		}

		cellSize = cellIdList->GetNumberOfIds();
		actualBound = GetBound(bifurcation,cellSize,this->Bound);
		localId = actualBound; // Start after the section included in the spine over the previous bifurcation.

		// Are there straight segments in the next iteration?
		if (actualBound < this->Bound || actualBound > this->Bound)
		{
			straightSegments = false;
		}
		else
		{
			straightSegments = true;
			if (cellSize-actualBound == this->Bound)
			{
				straightSegments = false;
			}
		}
	}

	output->SetLines(segments);

	return 1;

}


