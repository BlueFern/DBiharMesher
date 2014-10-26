/**
 * Program: vtkScalarRadiiToVectorsFilter.
 */

#include <vtkObjectFactory.h>
#include <vtkCallbackCommand.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>

#include "vtkScalarRadiiToVectorsFilter.h"

vtkStandardNewMacro(vtkScalarRadiiToVectorsFilter);

vtkScalarRadiiToVectorsFilter::vtkScalarRadiiToVectorsFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkScalarRadiiToVectorsFilter::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{

	// get the input and output
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	// input->GetPoints()->InsertNextPoint(1.0, 1.0, 1.0);

	output->ShallowCopy(input);

	// this->UpdateProgress(static_cast<double>(dim)/static_cast<double>(numDims));

	return 1;
}

void vtkScalarRadiiToVectorsFilter::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}

void vtkScalarRadiiToVectorsFilter::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkScalarRadiiToVectorsFilter* filter = static_cast<vtkScalarRadiiToVectorsFilter *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}
