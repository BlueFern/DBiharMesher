/**
 * This filter subdivides a quad mesh into smaller quads (roughly) of a specified size. The average
 * length and height of the existing quads are used to calculate how many new cells will fit inside
 * each one, which are then passed to vtkSubdivideQuadFilter.
 * Returned is a complete mesh with new quads, as vtkPolyData.
 */

#include <algorithm>
#include <cmath>
#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkCallbackCommand.h>
#include <vtkMath.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkAppendPolyData.h>

#include "vtkSubdivideMeshDynamic.h"
#include "vtkSubdivideQuadFilter.h"

vtkStandardNewMacro(vtkSubdivideMeshDynamic);

vtkSubdivideMeshDynamic::vtkSubdivideMeshDynamic()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Height = 0.0;
	this->Length = 0.0;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkSubdivideMeshDynamic::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	if (this->Height == 0.0 || this->Length == 0.0)
	{
		vtkErrorMacro("Must set both Columns and Rows to positive numbers.");
		exit(EXIT_FAILURE);
	}

	vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
	vtkSmartPointer<vtkIdList> pointsList = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	pointsList->SetNumberOfIds(4);
	points->SetNumberOfPoints(4);
	double p0[3], p1[3], p2[3], p3[3] = {0.0};
	int rows, columns = 0;
	double distance0, distance1, distance2, distance3 = 0.0;
	double averageDistance0, averageDistance1 = 0.0;
	int numberOfQuads = input->GetNumberOfCells();

	// Loop over every quad in input
	for (int i = 0; i < numberOfQuads; i++)
	{
		vtkSmartPointer<vtkSubdivideQuadFilter> subdivideQuadFilter = vtkSmartPointer<vtkSubdivideQuadFilter>::New();
		vtkSmartPointer<vtkPolyData> pointsData = vtkSmartPointer<vtkPolyData>::New();
		input->GetCellPoints(i,pointsList);
		points->SetPoint(0,input->GetPoint(pointsList->GetId(0)));
		points->SetPoint(1,input->GetPoint(pointsList->GetId(1)));
		points->SetPoint(2,input->GetPoint(pointsList->GetId(2)));
		points->SetPoint(3,input->GetPoint(pointsList->GetId(3)));

		points->GetPoint(0, p0);
		points->GetPoint(1, p1);
		points->GetPoint(2, p2);
		points->GetPoint(3, p3);

		distance0 = sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
		distance1 = sqrt(vtkMath::Distance2BetweenPoints(p3, p2));
		averageDistance0 = (distance0 + distance1) / 2.0;
		distance2 = sqrt(vtkMath::Distance2BetweenPoints(p0, p3));
		distance3 = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
		averageDistance1 = (distance2 + distance3) / 2.0;

		columns = round(averageDistance1 / this->Length);
		rows = round(averageDistance0 / this->Height);

		if (rows == 0)
		{
			rows++;
		}
		if (columns == 0)
		{
			columns++;
		}

		pointsData->SetPoints(points);
		subdivideQuadFilter->SetInputData(pointsData);

		subdivideQuadFilter->SetColumns(columns);
		subdivideQuadFilter->SetRows(rows);

		subdivideQuadFilter->Update();
		appendPolyData->AddInputData(subdivideQuadFilter->GetOutput());
	}
	appendPolyData->Update();
	output->ShallowCopy(appendPolyData->GetOutput());
	return 1;
}

void vtkSubdivideMeshDynamic::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Height: " << this->Height << " mm\n";
	os << indent << "Length: " << this->Length << " mm\n";
}

void vtkSubdivideMeshDynamic::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkSubdivideMeshDynamic* filter = static_cast<vtkSubdivideMeshDynamic *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}
