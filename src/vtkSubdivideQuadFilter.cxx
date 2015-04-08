#include <algorithm>

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkCallbackCommand.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkAppendPoints.h>
#include <vtkAppendPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>

#include "vtkSubdivideQuadFilter.h"

vtkStandardNewMacro(vtkSubdivideQuadFilter);

vtkSubdivideQuadFilter::vtkSubdivideQuadFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Columns = 0;
	this->Rows = 0;

	this->ShowProgress = false;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

/**
 * Copy the values from array1 into array2.
 */
void vtkSubdivideQuadFilter::copyPointsArray(double* array1, double* array2)
{
	for (int i = 0; i < 3; i++)
	{
		array2[i] = array1[i];
	}
}

int vtkSubdivideQuadFilter::RequestData(vtkInformation *vtkNotUsed(request),
		vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	// Test input data + params.

	if (this->Columns == 0 || this->Rows == 0)
	{
		vtkErrorMacro("Must set both Columns and Rows to positive integers.");
		exit(EXIT_FAILURE);
	}


	vtkSmartPointer<vtkPoints> topEdge = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> bottomEdge = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPoints> splineInputPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkAppendPoints> appendPoints = vtkSmartPointer<vtkAppendPoints>::New();
	vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();

	splineInputPoints->SetNumberOfPoints(2);

	if (this->Columns == 1)
	{
		topEdge->InsertPoint(0, input->GetPoint(0));
		topEdge->InsertPoint(1, input->GetPoint(3));

		bottomEdge->InsertPoint(0, input->GetPoint(1));
		bottomEdge->InsertPoint(1, input->GetPoint(2));
	}
	else
	{
		splineInputPoints->SetPoint(0, input->GetPoint(0));
		splineInputPoints->SetPoint(1, input->GetPoint(3));

		vtkSmartPointer<vtkParametricSpline> parametricSpline = vtkSmartPointer<vtkParametricSpline>::New();

		parametricSpline->SetPoints(splineInputPoints);
		// Setting these constraints ensures regular interval sampling of output points.
		parametricSpline->SetLeftConstraint(2);
		parametricSpline->SetRightConstraint(2);

		vtkSmartPointer<vtkParametricFunctionSource> splinePointsSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
		splinePointsSource->SetParametricFunction(parametricSpline);
		splinePointsSource->SetUResolution(this->Columns);
		splinePointsSource->Update();

		topEdge = splinePointsSource->GetOutput()->GetPoints();

		splineInputPoints->SetPoint(0,input->GetPoint(1));
		splineInputPoints->SetPoint(1,input->GetPoint(2));

		vtkSmartPointer<vtkParametricSpline> parametricSpline2 = vtkSmartPointer<vtkParametricSpline>::New();
		parametricSpline2->SetPoints(splineInputPoints);
		// Setting these constraints ensures regular interval sampling of output points.
		parametricSpline2->SetLeftConstraint(2);
		parametricSpline2->SetRightConstraint(2);

		vtkSmartPointer<vtkParametricFunctionSource> splinePointsSource2 = vtkSmartPointer<vtkParametricFunctionSource>::New();
		splinePointsSource2->SetParametricFunction(parametricSpline2);
		splinePointsSource2->SetUResolution(this->Columns);
		splinePointsSource2->Update();

		bottomEdge = splinePointsSource2->GetOutput()->GetPoints();
	}

	this->UpdateProgress(static_cast<double>(1) / static_cast<double>(this->Columns + 1));
	int test = topEdge->GetNumberOfPoints();

	// Build points down columns.

	if (this->Rows > 1)
	{
		for (int i = 0; i < this->Columns + 1; i++)
		{

			splineInputPoints->SetPoint(0,topEdge->GetPoint(i));
			splineInputPoints->SetPoint(1,bottomEdge->GetPoint(i));

			vtkSmartPointer<vtkParametricSpline> parametricSplineCol = vtkSmartPointer<vtkParametricSpline>::New();
			parametricSplineCol->SetPoints(splineInputPoints);
			parametricSplineCol->SetLeftConstraint(2);
			parametricSplineCol->SetRightConstraint(2);

			vtkSmartPointer<vtkParametricFunctionSource> sps = vtkSmartPointer<vtkParametricFunctionSource>::New();
			sps->SetParametricFunction(parametricSplineCol);
			sps->SetUResolution(this->Rows);
			sps->Update();

			vtkSmartPointer<vtkPolyData> pointSet = vtkSmartPointer<vtkPolyData>::New();
			pointSet->SetPoints(sps->GetOutput()->GetPoints());
			appendPoints->AddInputData(pointSet);

			this->UpdateProgress(static_cast<double>(i + 2) / static_cast<double>(this->Columns + 1));

		}
	}
	else
	{
		// No new rows to create, just use top and bottom edges.
		vtkSmartPointer<vtkPoints> allPoints = vtkSmartPointer<vtkPoints>::New();

		// Necessary for the correct ordering of points.
		for (int i = 0; i < topEdge->GetNumberOfPoints(); i++)
		{
			allPoints->InsertNextPoint(topEdge->GetPoint(i));
			allPoints->InsertNextPoint(bottomEdge->GetPoint(i));

		}

		vtkSmartPointer<vtkPolyData> pointSet = vtkSmartPointer<vtkPolyData>::New();
		pointSet->SetPoints(allPoints);
		appendPoints->AddInputData(pointSet);

	}

	appendPoints->Update();

	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetPoints(appendPoints->GetOutput()->GetPoints());
	structuredGrid->SetDimensions(this->Rows + 1, this->Columns + 1, 1);

	vtkSmartPointer<vtkStructuredGridGeometryFilter> structuredGridGeomFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	structuredGridGeomFilter->SetInputData(structuredGrid);

	int extent[6] = {0, this->Rows + 1, 0, this->Columns + 1, 0, 0};
	structuredGridGeomFilter->SetExtent(extent);
	structuredGridGeomFilter->Update();

	output->ShallowCopy(structuredGridGeomFilter->GetOutput());
	return 1;
}

void vtkSubdivideQuadFilter::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkSubdivideQuadFilter* filter = static_cast<vtkSubdivideQuadFilter *>(caller);
	if(filter->ShowProgress)
	{
		cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
	}
}

void vtkSubdivideQuadFilter::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Number of rows: " << this->Rows << "\n";
	os << indent << "Number of columns: " << this->Columns << "\n";
}
