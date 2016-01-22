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

#include "vtkSubdivideQuadBrick.h"
#include "vtkDbiharStatic.h"

vtkStandardNewMacro(vtkSubdivideQuadBrick);

vtkSubdivideQuadBrick::vtkSubdivideQuadBrick()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Columns = 0;
	this->Rows = 0;
	this->CellType = -1;
	this->ShowProgress = false;
	this->Filled = false;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkSubdivideQuadBrick::RequestData(vtkInformation *vtkNotUsed(request),
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

	if (this->Rows == 1)
	{
		topEdge->InsertPoint(0, input->GetPoint(0));
		topEdge->InsertPoint(1, input->GetPoint(1));

		bottomEdge->InsertPoint(0, input->GetPoint(3));
		bottomEdge->InsertPoint(1, input->GetPoint(2));
	}
	else
	{
		this->Rows *= 2; // Extra points for brick tessellation.

		splineInputPoints->SetPoint(0, input->GetPoint(0));
		splineInputPoints->SetPoint(1, input->GetPoint(3));

		vtkSmartPointer<vtkParametricSpline> parametricSpline = vtkSmartPointer<vtkParametricSpline>::New();

		parametricSpline->SetPoints(splineInputPoints);
		// Setting these constraints ensures regular interval sampling of output points.
		parametricSpline->SetLeftConstraint(2);
		parametricSpline->SetRightConstraint(2);

		vtkSmartPointer<vtkParametricFunctionSource> splinePointsSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
		splinePointsSource->SetParametricFunction(parametricSpline);
		splinePointsSource->SetUResolution(this->Rows);
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
		splinePointsSource2->SetUResolution(this->Rows);
		splinePointsSource2->Update();

		bottomEdge = splinePointsSource2->GetOutput()->GetPoints();
	}

	this->UpdateProgress(static_cast<double>(1) / static_cast<double>(this->Rows + 1));

	// Build points down columns.

	if (this->Columns > 1)
	{
		for (int i = 0; i < this->Rows + 1; i++)
		{

			splineInputPoints->SetPoint(0,topEdge->GetPoint(i));
			splineInputPoints->SetPoint(1,bottomEdge->GetPoint(i));

			vtkSmartPointer<vtkParametricSpline> parametricSplineCol = vtkSmartPointer<vtkParametricSpline>::New();
			parametricSplineCol->SetPoints(splineInputPoints);
			parametricSplineCol->SetLeftConstraint(2);
			parametricSplineCol->SetRightConstraint(2);

			vtkSmartPointer<vtkParametricFunctionSource> sps = vtkSmartPointer<vtkParametricFunctionSource>::New();
			sps->SetParametricFunction(parametricSplineCol);
			sps->SetUResolution(this->Columns);
			sps->Update();

			vtkSmartPointer<vtkPolyData> pointSet = vtkSmartPointer<vtkPolyData>::New();
			pointSet->SetPoints(sps->GetOutput()->GetPoints());
			appendPoints->AddInputData(pointSet);

			this->UpdateProgress(static_cast<double>(i + 2) / static_cast<double>(this->Rows + 1));
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

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();

	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	int baseId = this->Columns;
	int id = this->Columns;
	for (int i = 0; i < this->Rows / 2; i++)
	{
		for (int j = 0; j < this->Columns; j++)
		{
			vtkSmartPointer<vtkIdList> newCell = vtkSmartPointer<vtkIdList>::New();

			if (j % 2 == 0)
			{
				id -= this->Columns;
			}
			else
			{
				id += (this->Columns + 2);
			}

			if (i + 1 == this->Rows / 2 && j % 2 != 0) // Last iteration, do not extend cells beyond quad.
			{
				// Requires 6 points so that the later call to ReplaceCell is successful.
				newCell->InsertNextId(id);
				newCell->InsertNextId(id + (this->Columns + 1));
				newCell->InsertNextId(id + (this->Columns + 1));
				newCell->InsertNextId(id + (this->Columns + 2));
				newCell->InsertNextId(id + (this->Columns + 2));
				newCell->InsertNextId(id + 1);

			}
			else
			{
				newCell->InsertNextId(id);
				newCell->InsertNextId(id + (this->Columns + 1) * 2);
				newCell->InsertNextId(id + (this->Columns + 1) * 2);
				newCell->InsertNextId(id + (this->Columns + 1) * 2 + 1);
				newCell->InsertNextId(id + (this->Columns + 1) * 2 + 1);
				newCell->InsertNextId(id + 1);
			}

			cells->InsertNextCell(newCell);
		}
		baseId = baseId + (this->Columns + 1) * 2;
		id = baseId;
	}

	pd->SetPoints(appendPoints->GetOutput()->GetPoints());
	pd->SetPolys(cells);


	output->ShallowCopy(pd);
	return 1;
}

void vtkSubdivideQuadBrick::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkSubdivideQuadBrick* filter = static_cast<vtkSubdivideQuadBrick *>(caller);
	if(filter->ShowProgress)
	{
		cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
	}
}

void vtkSubdivideQuadBrick::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Number of rows: " << this->Rows << "\n";
	os << indent << "Number of columns: " << this->Columns << "\n";
}
