#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkCardinalSpline.h>

#include "vtkCentrelineData.h"

vtkStandardNewMacro(vtkCentrelineData);

// TODO: This class probably is better named as vtkCentrelineResampler.
// TODO: This class probably is better implemented as a filter.

vtkCentrelineData::vtkCentrelineData()
{
	polyData = vtkSmartPointer<vtkPolyData>::New();
}

void vtkCentrelineData::SetCentrelineData(vtkPolyData *centrelineData)
{
	// For each segment (vtkPolyLine) in centrelineData:
	// 1. Get segment length.
	// 2. Calculate the parametric resolution of the spline.
	// 3. Resample the segment.
	// 4. Resample the radii.
	// 4. Save points, lines, radii.

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	centrelineData->GetLines()->InitTraversal();
	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

	while(centrelineData->GetLines()->GetNextCell(cellIds))
	{
		// 1.
		vtkSmartPointer<vtkPoints> splineInputPoints = vtkSmartPointer<vtkPoints>::New();

		double length = 0;

		double p0[3];
		double p1[3];

		centrelineData->GetPoint(cellIds->GetId(0), p0);
		for(unsigned int i = 0; i < cellIds->GetNumberOfIds(); i++)
		{
			centrelineData->GetPoint(cellIds->GetId(i), p1);
			splineInputPoints->InsertNextPoint(p1);

			length += sqrt(vtkMath::Distance2BetweenPoints(p0, p1));

			memcpy(p0, p1, sizeof(double_t) * 3);
		}

		// 2.
		int resolution = (int)((length * this->unitsConversionFactor) / (ECMultiple * ECLength));

		// 3.
		vtkSmartPointer<vtkParametricSpline> parametricSpline = vtkSmartPointer<vtkParametricSpline>::New();
		parametricSpline->SetPoints(splineInputPoints);

		vtkSmartPointer<vtkParametricFunctionSource> splinePointsSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
		splinePointsSource->SetParametricFunction(parametricSpline);
		splinePointsSource->SetUResolution(resolution);
		splinePointsSource->Update();

		vtkPoints *splineOutputPoints = splinePointsSource->GetOutput()->GetPoints();

		unsigned int numOutputPoints = splineOutputPoints->GetNumberOfPoints();

		// 4.
		vtkSmartPointer<vtkCardinalSpline> radiiSpline = vtkSmartPointer<vtkCardinalSpline>::New();
		radiiSpline->SetParametricRange(0.0, length);

		centrelineData->GetPoint(cellIds->GetId(0), p0);
		double distance = 0.0;
		vtkDataArray *splineInputRadii = centrelineData->GetPointData()->GetScalars();
		for(unsigned int i = 0; i < cellIds->GetNumberOfIds(); i++)
		{
			centrelineData->GetPoint(cellIds->GetId(i), p1);
			distance += sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
			radiiSpline->AddPoint(distance, splineInputRadii->GetTuple(cellIds->GetId(i))[0]);

			memcpy(p0, p1, sizeof(double_t) * 3);
		}

		// 5.
		vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
		for(unsigned int i = 0; i < numOutputPoints; ++i)
		{
			vtkIdType pId = points->InsertNextPoint(splineOutputPoints->GetPoint(i));

			polyLine->GetPointIds()->InsertNextId(pId);

			double t = i * (length / (double)(numOutputPoints - 1));
			double radius = radiiSpline->Evaluate(t);
			radii->InsertNextTuple(&radius);
		}

		lines->InsertNextCell(polyLine);
	}

	polyData->SetPoints(points);
	polyData->SetLines(lines);
	polyData->GetPointData()->SetScalars(radii);
}

vtkPolyData *vtkCentrelineData::GetOutput()
{
	return polyData;
}

void vtkCentrelineData::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}
