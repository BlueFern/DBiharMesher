#include <cmath>
#include <algorithm>

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkCallbackCommand.h>

#include "vtkDbiharStatic.h"
#include "vtkPointsToMeshFilter.h"

vtkStandardNewMacro(vtkPointsToMeshFilter);


vtkPointsToMeshFilter::vtkPointsToMeshFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	// WARNING: VTK developers style guide says all class variables need to be initialised,
	// but in this case it indicates memory problems. This issue needs to be tested further.
	// this->Dimensions = vtkSmartPointer<vtkUnsignedIntArray>::New(); // Error.

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkPointsToMeshFilter::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	// Test Dimensions is not NULL and has been initialised.
	if (Dimensions->GetNumberOfTuples() < 2)
	{
		vtkErrorMacro("Require 3 or more tuples in Dimensions array (Got " << Dimensions->GetNumberOfTuples() << ").");
		exit(EXIT_FAILURE);
	}

	// Calculate number of patches depending on the number of items in the dimensions array.
	// If numPatches is 1, straight line segment. Otherwise a bifurcation (3 or more for a higher branching number).
	int numPatches = Dimensions->GetNumberOfTuples() - 1;

	// Test the number of points in the input matches the dimensions.
	int numQuads = 0;
	int numPoints = 0;

	for (int i = 1; i < Dimensions->GetNumberOfTuples(); i++)
	{
		numQuads += 2 * (Dimensions->GetValue(i) + 1);
	}
	if (numPatches > 1)
	{
		numQuads -= numPatches; // Removes overlapping points between branches.
	}
	numPoints = (Dimensions->GetValue(0) + 1) * numQuads;

	if (input->GetNumberOfPoints() != numPoints)
	{
		vtkErrorMacro("Number of input points (" << input->GetNumberOfPoints() << ") does not match the number of points specified by dimensions (" << numPoints << ").");
		exit(EXIT_FAILURE);
	}

	vtkSmartPointer<vtkPolyData> result = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkIntArray> cellData = vtkSmartPointer<vtkIntArray>::New();
	vtkSmartPointer<vtkCellArray> quads = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkIdList> quad = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	double point[3] = {0,0,0};
	int start = 0;
	int end = 0;
	int halfLoop = (Dimensions->GetValue(0) + 1);
	int numPointsLoop = 2 * Dimensions->GetValue(0); // Only used for quads.

	// Initial values for trunk. Will be updated at the end of each iteration for the next branch.
	int branchStart = 0;
	int reversedStart = numPoints - 1;
	int reversedEnd = reversedStart - halfLoop;
	int quadPosition = 0;
	int cellDataId = 0;

	for (int branch = 1; branch <= numPatches; branch++) // Looping over the number of branches (once if straight segment).
	{
		for (int i = 0; i < Dimensions->GetValue(branch) + 1; i++) // For each ring in a given branch.
		{
			start = branchStart + i * (Dimensions->GetValue(0) + 1);
			end = start + (Dimensions->GetValue(0) + 1);

			for (int j = start; j < end; j++) // Bottom half of ring.
			{
				input->GetPoint(j, point);
				points->InsertNextPoint(point);

				// Build quads. The points of the quads are not yet in the points array, but we know where they will be based on
				// values in the input dimensions array.
				if (i < Dimensions->GetValue(branch) && j + 1 < end) // Don't make quads on last loop.
				{
					quad->InsertUniqueId(quadPosition);
					quad->InsertUniqueId(quadPosition + 1);
					quad->InsertUniqueId(quadPosition + numPointsLoop + 1);
					quad->InsertUniqueId(quadPosition + numPointsLoop);
					quads->InsertNextCell(quad);
					cellData->InsertNextValue(cellDataId);
					quad->Reset();
					quadPosition++;
				}
			}

			reversedEnd = reversedStart - halfLoop;

			for (int j = reversedStart; j > reversedEnd; j--) // Top half of the ring.
			{

				if (j < reversedStart && j > reversedEnd + 1) // Duplicate points.
				{
					input->GetPoint(j, point);
					points->InsertNextPoint(point);
				}

				// Creating top half quads.
				if (i < Dimensions->GetValue(branch) && j > reversedEnd + 2)
				{
					quad->InsertUniqueId(quadPosition);
					quad->InsertUniqueId(quadPosition + 1);
					quad->InsertUniqueId(quadPosition + numPointsLoop + 1);
					quad->InsertUniqueId(quadPosition + numPointsLoop);
					quads->InsertNextCell(quad);
					cellData->InsertNextValue(cellDataId);
					quad->Reset();
					quadPosition++;
				}
			}

			// Connecting the last quad in the ring to the first.
			if (i < Dimensions->GetValue(branch)) // Again, no quads on the last ring.
			{
				quad->InsertUniqueId(quadPosition);
				quad->InsertUniqueId(quadPosition - (numPointsLoop - 1));
				quad->InsertUniqueId(quadPosition + 1);
				quad->InsertUniqueId(quadPosition + numPointsLoop);
				quads->InsertNextCell(quad);
				cellData->InsertNextValue(cellDataId);
				quad->Reset();
				quadPosition++;
			}
			reversedStart -= halfLoop;
		}

		// End early if in last iteration of loop.
		if (branch + 1 > numPatches)
		{
			break;
		}

		cellDataId++;

		// Find the starting points of the rings and quads in the next branch.

		branchStart = halfLoop * (Dimensions->GetValue(1) + 1);
		for (int k = 2; k <= branch; k++)
		{
			branchStart +=  2 * halfLoop * (Dimensions->GetValue(k) + 1);
		}

		// Duplicating points between branches is intended. Move back half a loop for every branch we've passed.
		branchStart -= branch * halfLoop;

		reversedStart = branchStart + (2 * halfLoop * (Dimensions->GetValue(branch + 1) + 1)) - 1;

		quadPosition += numPointsLoop;

	}
	cellData->SetName(vtkDbiharStatic::CELL_DATA_ARR_NAME);
	result->GetCellData()->SetScalars(cellData);
	result->SetPoints(points);
	result->SetPolys(quads);

	output->ShallowCopy(result);

	// Required to return 1 by VTK API.
	return 1;
}

void vtkPointsToMeshFilter::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	// this->Dimensions->PrintSelf(os, indent);
}

void vtkPointsToMeshFilter::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkPointsToMeshFilter* filter = static_cast<vtkPointsToMeshFilter *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}
