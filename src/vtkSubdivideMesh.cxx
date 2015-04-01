//#include <algorithm>
//#include <cmath>

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkCallbackCommand.h>


#include <vtkAppendPolyData.h>

#include "vtkSubdivideMesh.h"
#include "vtkSubdivideQuadFilter.h"

vtkStandardNewMacro(vtkSubdivideMesh);

vtkSubdivideMesh::vtkSubdivideMesh()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Columns = 0;
	this->Rows = 0;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkSubdivideMesh::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	if (this->Columns == 0.0 || this->Rows == 0.0)
	{
		vtkErrorMacro("Must set both Columns and Rows to positive numbers.");
		exit(EXIT_FAILURE);
	}

	int numberOfQuads = input->GetNumberOfCells();

	int percentOfTotal  = numberOfQuads / 10; // Every 10%, for progress function.
	int stage = 1;

	vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
	vtkSmartPointer<vtkIdList> pointsList = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	pointsList->SetNumberOfIds(4);
	points->SetNumberOfPoints(4);

	for (int i = 0; i < input->GetNumberOfCells(); i++)
	{
		vtkSmartPointer<vtkSubdivideQuadFilter> subdivideQuadFilter = vtkSmartPointer<vtkSubdivideQuadFilter>::New();
		vtkSmartPointer<vtkPolyData> pointsData = vtkSmartPointer<vtkPolyData>::New();
		input->GetCellPoints(i,pointsList);
		points->SetPoint(0,input->GetPoint(pointsList->GetId(0)));
		points->SetPoint(1,input->GetPoint(pointsList->GetId(1)));
		points->SetPoint(2,input->GetPoint(pointsList->GetId(2)));
		points->SetPoint(3,input->GetPoint(pointsList->GetId(3)));
		pointsData->SetPoints(points);
		subdivideQuadFilter->SetInputData(pointsData);
		subdivideQuadFilter->SetColumns(this->Columns);
		subdivideQuadFilter->SetRows(this->Rows);
		subdivideQuadFilter->Update();
		appendPolyData->AddInputData(subdivideQuadFilter->GetOutput());

		// Basic progress reporting.
		if ((i + 1) % percentOfTotal == 0)
		{
			this->UpdateProgress(static_cast<double>(stage++) / static_cast<double>(11)); // 10% + 1 as we start with stage at 1.
		}

	}
	appendPolyData->Update();

	output->ShallowCopy(appendPolyData->GetOutput());
	return 1;
}

void vtkSubdivideMesh::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Columns: " << this->Columns << "\n";
	os << indent << "Rows: " << this->Rows << "\n";
}

void vtkSubdivideMesh::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkSubdivideMesh* filter = static_cast<vtkSubdivideMesh *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}
