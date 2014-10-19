#include <map>

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
	// 0. Find shared point ids for each segment.
	// 1. Get segment length.
	// 2. Calculate the parametric resolution of the spline.
	// 3. Resample the segment.
	// 4. Resample the radii.
	// 4. Save points, lines, radii.

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

	// 0.

	// Mapping from branch to parent segments.
	std::map<vtkIdType, vtkIdType> sharedIds;

	unsigned int segmentId = 0;
	centrelineData->GetLines()->InitTraversal();
	while(centrelineData->GetLines()->GetNextCell(cellIds))
	{
		vtkSmartPointer<vtkIdList> lastId = vtkSmartPointer<vtkIdList>::New();
		lastId->InsertNextId(cellIds->GetId(cellIds->GetNumberOfIds() - 1));

		// WARNING: The code in the GetCellNeighbors method changes the traversal location in the vtkCellArray. In my view it is utterly retarded.
		vtkIdType traverseLocation = centrelineData->GetLines()->GetTraversalLocation();

		vtkSmartPointer<vtkIdList> neighbourCellIds = vtkSmartPointer<vtkIdList>::New();
		centrelineData->GetCellNeighbors(segmentId, lastId, neighbourCellIds);

		// Restore traversal location.
		centrelineData->GetLines()->SetTraversalLocation(traverseLocation);

		// Save the segment id in the shared ids map with the corresponding parent segment id.
		for(vtkIdType i = 0; i < neighbourCellIds->GetNumberOfIds(); i++)
		{
			sharedIds[neighbourCellIds->GetId(i)] = segmentId;
		}
		segmentId++;
	}

	// Mapping from segment to last point ids.
	std::map<vtkIdType, vtkIdType> lastPointIds;

	segmentId = 0;
	centrelineData->GetLines()->InitTraversal();
	while(centrelineData->GetLines()->GetNextCell(cellIds))
	{
		std::cout << centrelineData->GetLines()->GetTraversalLocation() << std::endl;

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
		for(vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
		{
			centrelineData->GetPoint(cellIds->GetId(i), p1);
			distance += sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
			radiiSpline->AddPoint(distance, splineInputRadii->GetTuple(cellIds->GetId(i))[0]);

			memcpy(p0, p1, sizeof(double_t) * 3);
		}

		// 5.

		// TODO: Review this: somehow the distances between points end up being slightly irregular even after resampling.
		// splineOutputPoints->GetPoint(0, p0);

		vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
		for(unsigned int i = 0; i < numOutputPoints; ++i)
		{
			// splineOutputPoints->GetPoint(i, p1);
			// std::cout << sqrt(vtkMath::Distance2BetweenPoints(p0, p1)) << std::endl;
			// memcpy(p0, p1, sizeof(double_t) * 3);

			vtkIdType pId = -1;

			// Check if this branch has a parent segment.
			if(i == 0)
			{
				// See if parent for this branch exists.
				std::map<vtkIdType, vtkIdType>::iterator pos1 = sharedIds.find(segmentId);
				if(pos1 != sharedIds.end())
				{
					vtkIdType lineId = pos1->second;

					// Find the last inserted point of the parent segment.
					std::map<vtkIdType, vtkIdType>::iterator pos2 = lastPointIds.find(lineId);
					if(pos2 != lastPointIds.end())
					{
						pId = pos2->second;
					}
					else
					{
						// This should not normally occur.
						vtkErrorWithObjectMacro(this, "Something is rotten in the state of Denmark: " << __FILE__ << ":" << __LINE__);
					}
				}
			}
			// No shared point ids found.
			if(pId == -1)
			{
				pId = points->InsertNextPoint(splineOutputPoints->GetPoint(i));
			}

			polyLine->GetPointIds()->InsertNextId(pId);

			double t = i * (length / (double)(numOutputPoints - 1));
			double radius = radiiSpline->Evaluate(t);
			radii->InsertNextTuple(&radius);

			// Remember the last stored id for this segment.
			if(i == numOutputPoints - 1)
			{
				lastPointIds[segmentId] = pId;
			}
		}

		lines->InsertNextCell(polyLine);
		segmentId++;
	}

#if 0
	for (std::map<vtkIdType,vtkIdType>::iterator it=sharedIds.begin(); it!=sharedIds.end(); ++it)
		std::cout << it->first << " => " << it->second << '\n';

	for (std::map<vtkIdType,vtkIdType>::iterator it=lastPointIds.begin(); it!=lastPointIds.end(); ++it)
		std::cout << it->first << " => " << it->second << '\n';
#endif

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
