#include <map>

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkPolyLine.h>
#include <vtkCallbackCommand.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkCardinalSpline.h>

#include "vtkDbiharStatic.h"
#include "vtkCentrelineResampler.h"

vtkStandardNewMacro(vtkCentrelineResampler);

vtkCentrelineResampler::vtkCentrelineResampler()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->EdgeLength = 0;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkCentrelineResampler::RequestData(vtkInformation *vtkNotUsed(request),
										vtkInformationVector **inputVector,
										vtkInformationVector *outputVector)
{
	vtkPolyData *input = vtkPolyData::GetData(inputVector[0],0);
	input->BuildLinks();

	vtkPolyData *output = vtkPolyData::GetData(outputVector,0);

	// For each segment (vtkPolyLine) in centrelineData:
	// 0. Find shared point ids for each segment.
	// 1. Get segment length.
	// 2. Calculate the parametric resolution of the spline.
	// 3. Resample the segment.
	// 4. Resample the radii.
	// 5. Save points, lines, radii.

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkDoubleArray> radii = vtkSmartPointer<vtkDoubleArray>::New();
	radii->SetName(vtkDbiharStatic::RADII_SCALARS_ARR_NAME);

	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

	// 0.

	// Mapping from branch to parent segments.
	std::map<vtkIdType, vtkIdType> sharedIds;

	unsigned int segmentId = 0;
	input->GetLines()->InitTraversal();
	while(input->GetLines()->GetNextCell(cellIds))
	{
		vtkSmartPointer<vtkIdList> lastId = vtkSmartPointer<vtkIdList>::New();
		lastId->InsertNextId(cellIds->GetId(cellIds->GetNumberOfIds() - 1));

		// WARNING: The code in the GetCellNeighbors method changes the traversal location in the vtkCellArray. In my view it is utterly retarded.
		vtkIdType traverseLocation = input->GetLines()->GetTraversalLocation();

		vtkSmartPointer<vtkIdList> neighbourCellIds = vtkSmartPointer<vtkIdList>::New();
		input->GetCellNeighbors(segmentId, lastId, neighbourCellIds);

		// Restore traversal location.
		input->GetLines()->SetTraversalLocation(traverseLocation);

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
	input->GetLines()->InitTraversal();
	int numLines = input->GetLines()->GetNumberOfCells() + 1; // For progress function.
	int counter = 1;
	while(input->GetLines()->GetNextCell(cellIds))
	{
		// 1.
		vtkSmartPointer<vtkPoints> splineInputPoints = vtkSmartPointer<vtkPoints>::New();

		double length = 0;

		double p0[3];
		double p1[3];

		input->GetPoint(cellIds->GetId(0), p0);
		for(unsigned int i = 0; i < cellIds->GetNumberOfIds(); i++)
		{
			input->GetPoint(cellIds->GetId(i), p1);
			splineInputPoints->InsertNextPoint(p1);

			length += sqrt(vtkMath::Distance2BetweenPoints(p0, p1));

			memcpy(p0, p1, sizeof(double_t) * 3);
		}


		// 2. Make sure each cell contains odd number of points.
		int resolution = (int)(length / this->EdgeLength);
		// If the initial value for the resolution is even, adjust it by 1 in the direction of the least error.
		if(resolution & 1)
		{
			double tmp = std::fabs(length - resolution * this->EdgeLength);
			// If the remainder is less then half of the given edge length, decrement the resolution value. Otherwise increment it.
			if(tmp < this->EdgeLength / 2.0)
			{
				resolution--;
			}
			else
			{
				resolution++;
			}
		}
		// TODO: Print the error of the actual edge length using vtkWarningMacro.

		// 3.
		vtkSmartPointer<vtkParametricSpline> parametricSpline = vtkSmartPointer<vtkParametricSpline>::New();
		parametricSpline->SetPoints(splineInputPoints);
		// Setting these constraints ensures regular interval sampling of output points.
		parametricSpline->SetLeftConstraint(2);
		parametricSpline->SetRightConstraint(2);

		vtkSmartPointer<vtkParametricFunctionSource> splinePointsSource = vtkSmartPointer<vtkParametricFunctionSource>::New();
		splinePointsSource->SetParametricFunction(parametricSpline);
		splinePointsSource->SetUResolution(resolution);
		splinePointsSource->Update();

		vtkPoints *splineOutputPoints = splinePointsSource->GetOutput()->GetPoints();

		unsigned int numOutputPoints = splineOutputPoints->GetNumberOfPoints();

		// 4.
		vtkSmartPointer<vtkCardinalSpline> radiiSpline = vtkSmartPointer<vtkCardinalSpline>::New();
		radiiSpline->SetParametricRange(0.0, length);

		input->GetPoint(cellIds->GetId(0), p0);
		double distance = 0.0;
		vtkDataArray *splineInputRadii = input->GetPointData()->GetScalars();
		for(vtkIdType i = 0; i < cellIds->GetNumberOfIds(); i++)
		{
			input->GetPoint(cellIds->GetId(i), p1);
			distance += sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
			radiiSpline->AddPoint(distance, splineInputRadii->GetTuple(cellIds->GetId(i))[0]);

			memcpy(p0, p1, sizeof(double_t) * 3);
		}
		;
		// 5.
		vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
		for(unsigned int i = 0; i < numOutputPoints; ++i)
		{
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
						vtkErrorMacro("Something is rotten in the state of Denmark: " << __FILE__ << ":" << __LINE__);
						exit(EXIT_FAILURE);
					}
				}
			}
			// No shared point ids found.
			if(pId == -1)
			{
				pId = points->InsertNextPoint(splineOutputPoints->GetPoint(i));

				double t = i * (length / (double)(numOutputPoints - 1));
				double radius = radiiSpline->Evaluate(t);
				radii->InsertNextTuple(&radius);
			}

			polyLine->GetPointIds()->InsertNextId(pId);

			// Remember the last stored id for this segment.
			if (i == numOutputPoints - 1)
			{
				lastPointIds[segmentId] = pId;
			}
		}

		if(!polyLine->GetNumberOfPoints() & 1)
		{
			vtkErrorMacro("Something is rotten in the state of Denmark:" << __FILE__ << ":" << __LINE__);
			exit(EXIT_FAILURE);
		}

		lines->InsertNextCell(polyLine);
		segmentId++;
		this->UpdateProgress(static_cast<double>(counter++) / static_cast<double>(numLines));
	}

#if 0
	for (std::map<vtkIdType,vtkIdType>::iterator it=sharedIds.begin(); it!=sharedIds.end(); ++it)
		std::cout << it->first << " => " << it->second << '\n';

	for (std::map<vtkIdType,vtkIdType>::iterator it=lastPointIds.begin(); it!=lastPointIds.end(); ++it)
		std::cout << it->first << " => " << it->second << '\n';
#endif

	output->SetPoints(points);
	output->SetLines(lines);
	output->GetPointData()->SetScalars(radii);

	return 1;
}


void vtkCentrelineResampler::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "\nEdge Length: " << this->EdgeLength << "\n";
}

void vtkCentrelineResampler::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkCentrelineResampler* filter = static_cast<vtkCentrelineResampler *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}
