/*
 * Program: vtkPointsToMeshFilter.
 */

#include <cmath>
#include <algorithm>

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>

#include <vtkCallbackCommand.h>

#include "vtkPointsToMeshFilter.h"

vtkStandardNewMacro(vtkPointsToMeshFilter);

vtkPointsToMeshFilter::vtkPointsToMeshFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	// WARNING: VTK developers style guide says all class variables need to be initialised,
	// but in this case it indicates memory problems. This issue needs to be tested further.
	// this->Dimensions = vtkSmartPointer<vtkUnsignedIntArray>::New();

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

	// Test the number of items is 3 (straight segment) or 7 (bifurcation). More cases?

	// Calculate number of patches depending on the number of items in the dimensions array.
	int numPatches = (Dimensions->GetNumberOfTuples() - 1) / 2;

	// Decide here if dealing with bifurcation depending on the number of patches.

	// Test the number of points in the input matches the dimensions.
	int numQuads = 0;
	for(int i = 1; i < Dimensions->GetNumberOfTuples(); i++)
	{
		numQuads += Dimensions->GetValue(i);
	}

	int numPoints = (Dimensions->GetValue(0) + 1) * (numQuads + numPatches);

	if(input->GetNumberOfPoints() != numPoints)
	{
		vtkErrorMacro("Number of input points (" << input->GetNumberOfPoints() << ") does not match the number of points specified by dimensions (" << numPoints << ").");
		exit(EXIT_FAILURE);
	}

	// Produce polydata mesh.

	output->ShallowCopy(input);

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
