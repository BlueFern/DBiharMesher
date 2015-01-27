
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
#include <vtkCell.h>

#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkAppendPoints.h>
#include <vtkAppendPolyData.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkSubdivideQuadFilter.h>


vtkStandardNewMacro(vtkSubdivideQuadFilter);


vtkSubdivideQuadFilter::vtkSubdivideQuadFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Columns = 0;
	this->Rows = 0;
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

	vtkSmartPointer<vtkPoints> splineInputPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkAppendPoints> appendPoints = vtkSmartPointer<vtkAppendPoints>::New();
	vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();

	splineInputPoints->SetNumberOfPoints(2);

	splineInputPoints->SetPoint(0,input->GetPoint(0));
	splineInputPoints->SetPoint(1,input->GetPoint(3));

	vtkSmartPointer<vtkParametricSpline> parametricSpline = vtkSmartPointer<vtkParametricSpline>::New();

	parametricSpline->SetPoints(splineInputPoints);
	// Setting these constraints ensures regular interval sampling of output points.
	parametricSpline->SetLeftConstraint(2);
	parametricSpline->SetRightConstraint(2);

	vtkSmartPointer<vtkParametricFunctionSource> splinePointsSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
	splinePointsSource->SetParametricFunction(parametricSpline);
	splinePointsSource->SetUResolution(this->Columns);
	splinePointsSource->Update();

	vtkSmartPointer<vtkPoints> topEdge = splinePointsSource->GetOutput()->GetPoints();

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

	vtkSmartPointer<vtkPoints> bottomEdge = splinePointsSource2->GetOutput()->GetPoints();

	// Build points down columns.
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

void vtkSubdivideQuadFilter::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Number of rows: " << this->Rows << "\n";
	os << indent << "Number of columns: " << this->Columns << "\n";
}
