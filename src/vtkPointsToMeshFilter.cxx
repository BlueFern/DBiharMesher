/**
 * This filter builds a quadrilateral mesh using points specified in input vtkPolyData. This input
 * contains no cells or attribute data, only points that form a 3D model of a straight segment
 * or bifurcation. The order of the points is important; they follow the spines around
 * a bifurcation or up then down a straight segment in half loops/circles.
 *
 * The Dimensions array gives details about the particular input vtkPolyData. The first
 * element specifies the number of interior quads in each half ring (so is one less than
 * the number of points in it). The next one or more elements are the number of rows of
 * quads that will fit in a branch (so one less than the number of half rings in it).
 * (There would therefore be Dimensions->GetValue(0) * Dimensions->GetValue(b) quads on
 * branch 'b'.)
 *
 * vtkPolyData is returned with the quads as cells and most of the points from the input.
 * Duplicate points between connecting halves are omitted, but the points between connected
 * branches are included for each branch it is associated with. This association is actually
 * fictitious, as there are no cells or point data, but is nonetheless made for simplicity.
 * Keeping track of what branch (iteration of main loop) added a particular point created
 * unnecessary difficulty when building the quads.
 * The order of these points is also different. Instead of half rings the points move entirely
 * around the model. If we're not dealing with straight segments, points start at the bifurcation
 * and move out to the end of the branch.
 *
 * A vtkIntArray is created and used as cell data for the quads. It serves to
 * reference what trunk a particular cell belongs to; so for each quad created a value
 * is put into the array (0 if the cell belongs to the first trunk, 1 for the next, etc.).
 *
 */

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

#include "vtkPointsToMeshFilter.h"

vtkStandardNewMacro(vtkPointsToMeshFilter);

const char *vtkPointsToMeshFilter::CELL_DATA_ARR_NAME = {"branchId"};

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

		reversedStart = branchStart + (2 * halfLoop * (Dimensions->GetValue(branch+1) + 1)) - 1;

		quadPosition += numPointsLoop;

	}
	cellData->SetName(vtkPointsToMeshFilter::CELL_DATA_ARR_NAME);
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
