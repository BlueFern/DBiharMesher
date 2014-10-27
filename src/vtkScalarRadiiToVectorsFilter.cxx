/**
 * Program: vtkScalarRadiiToVectorsFilter.
 */

#include <map>

#include <vtkObjectFactory.h>
#include <vtkCallbackCommand.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>

#include <vtkMath.h>
#include <vtkTransform.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>

#include "vtkScalarRadiiToVectorsFilter.h"

#define PRINT_DEBUG 0

vtkStandardNewMacro(vtkScalarRadiiToVectorsFilter);

const char *vtkScalarRadiiToVectorsFilter::RADII_ARR_NAME = {"radiiVectors"};

void DoubleCross(const double v0[3], const double c0[3], const double v1[3], double c1[3])
{
	vtkMath::Cross(c0, v0, c1);
	vtkMath::Cross(v1, c1, c1);
}

// This is taken from vtkMath class in VTK Nightly in October 2014.
double AngleBetweenVectors(const double v1[3], const double v2[3])
{
  double cross[3];
  vtkMath::Cross(v1, v2, cross);
  return atan2(vtkMath::Norm(cross), vtkMath::Dot(v1, v2));
}

vtkScalarRadiiToVectorsFilter::vtkScalarRadiiToVectorsFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	angleTolerance = 1.0e-5;

	// TODO: Implement progress updates in RequestData.
	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkScalarRadiiToVectorsFilter::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	vtkSmartPointer<vtkCellArray> lines = input->GetLines();
	lines->InitTraversal();

	// Obtain tree structure (bifurcations and terminals) from the centreline tree.
	std::map<vtkIdType, std::vector<vtkIdType> > treeInfo;

	// Find treeInfo and terminals.
	for(vtkIdType lineId = 0; lineId < lines->GetNumberOfCells(); lineId++)
	{
		vtkSmartPointer<vtkIdList> lineIds = vtkSmartPointer<vtkIdList>::New();
		lines->GetNextCell(lineIds);
		vtkSmartPointer<vtkIdList> lastId = vtkSmartPointer<vtkIdList>::New();
		lastId->InsertNextId(lineIds->GetId(lineIds->GetNumberOfIds() - 1));

		// WARNING: The code in the GetCellNeighbors method changes the traversal location in the vtkCellArray. In my view it is utterly retarded.
		vtkIdType traverseLocation = input->GetLines()->GetTraversalLocation();

		vtkSmartPointer<vtkIdList> neighbourCellIds = vtkSmartPointer<vtkIdList>::New();
		input->GetCellNeighbors(lineId, lastId, neighbourCellIds);

		// Restore traversal location.
		input->GetLines()->SetTraversalLocation(traverseLocation);

		std::cout << "Line " << lineId << " shares last point with " << neighbourCellIds->GetNumberOfIds() << " cells." << std::endl;

		std::vector<vtkIdType> idList;

		// The assumption is that non-bifurcating segments are not divided in segments.
		if(neighbourCellIds->GetNumberOfIds() == 1)
		{
			vtkWarningMacro("Found a point shared between only two cells, which breaks the assumption of non-bifurcating segments integrity.");
		}

		for(vtkIdType pId = 0; pId < neighbourCellIds->GetNumberOfIds(); pId++)
		{
			idList.push_back(neighbourCellIds->GetId(pId));
		}
		treeInfo[lineId] = idList;
	}

#if PRINT_DEBUG
	// Print the map.
	for(std::map<vtkIdType, std::vector<vtkIdType> >::iterator it = treeInfo.begin(); it != treeInfo.end(); ++it)
	{
		std::vector<vtkIdType> v = it->second;
		std::cout << it->first << " => ";
		for(std::vector<vtkIdType>::iterator i = v.begin(); i != v.end(); ++i)
		{
			std::cout << *i << " ";
		}
		std::cout << endl;
	}
#endif

	// Array 0 is for upward radii.
	vtkSmartPointer<vtkDoubleArray> radiiVectors = vtkSmartPointer<vtkDoubleArray>::New();
	radiiVectors->SetNumberOfComponents(3);
	radiiVectors->SetNumberOfTuples(input->GetNumberOfPoints());
	radiiVectors->SetName(RADII_ARR_NAME);
	for(int i = 0; i < radiiVectors->GetNumberOfTuples(); i++)
	{
		radiiVectors->SetTuple3(i, 0, 0, 0);
	}

	double zero[3] = {0};
	double v0[3];
	double v1[3];
	double c0[3];
	double c1[3];

	// Keeping track of the average vectors at bifurcations.
	// This is a somewhat clunky way of doing so, which a map from
	// bifurcation point ids to indices in the bifurcation vectors array.
	std::map<vtkIdType, vtkIdType> lineIdToAverageVectorMap;
	vtkSmartPointer<vtkDoubleArray> bifurcationVectors = vtkSmartPointer<vtkDoubleArray>::New();
	bifurcationVectors->SetNumberOfComponents(3);

	// Process the bifurcations.
	for(std::map<vtkIdType, std::vector<vtkIdType> >::iterator it = treeInfo.begin(); it != treeInfo.end(); ++it)
	{
		// Skip all terminals
		if(it->second.size() == 0)
		{
			continue;
		}

		// Find radius vector c0 at this bifurcation.
		GetDirectionVector(it->second[0], d_start, v0);
		vtkMath::MultiplyScalar(v0, -1);
		GetDirectionVector(it->second[1], d_start, v1);
		vtkMath::MultiplyScalar(v1, -1);

		// c0 orthogonal to v0 and v1.
		vtkMath::Cross(v0, v1, c0);
		vtkMath::Normalize(c0);

		// Global point id.
		vtkIdType bifId = GetGlobalPointId(it->first, d_end);

		// Store bifurcation vector c0.
		radiiVectors->SetTuple(bifId, c0);

		// v0 here is an average of the two vectors.
		vtkMath::Add(v0, v1, v0);
		vtkMath::MultiplyScalar(v0, 0.5);

		// Store the average vector at the start of this bifurcation.
		lineIdToAverageVectorMap[it->second[0]] = bifurcationVectors->InsertNextTuple(v0);
		lineIdToAverageVectorMap[it->second[1]] = bifurcationVectors->InsertNextTuple(v0);

		// Rotation angle for radius at every point, to make sure there is a
		// gradual twist between start and end radius vector for each segment.
		double twistAngle = 0.0;
		double twistDirection = 1;

		// Figure out the twist angle between bifurcations.
		// Skip the first bifurcation.
		if(it->first > treeInfo.begin()->first)
		{
			vtkIdType parentBifId = GetGlobalPointId(it->first, d_start);
			bifurcationVectors->GetTuple((lineIdToAverageVectorMap[it->first]), v1);

			DoubleCross(v0, c0, v1, c1);
			vtkMath::Normalize(c1);

			double tmp[3];
			// c0 from parent bifurcation.
			radiiVectors->GetTuple(parentBifId, tmp);

			// Figure out the angle between this bifurcation, and the closest upstearm bifurcation.
			twistAngle = vtkMath::DegreesFromRadians(AngleBetweenVectors(c1, tmp));

			// Figure out twist direction.
			// Probably can use double cross instead of a transform.
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			transform->RotateWXYZ(twistAngle, v1);
			transform->TransformPoint(tmp, tmp);

			if(vtkMath::DegreesFromRadians(AngleBetweenVectors(tmp, c1)) < angleTolerance)
			{
				twistDirection = -1;
			}
		}

		// Process the rest of the line.
		vtkSmartPointer<vtkIdList> lineIds = GetLineIds(it->first);

		// Insert radii for the rest of the line.
		// Traversing the centreline in reverse direction.
		for(vtkIdType pId = lineIds->GetNumberOfIds() - 2; pId >= 0; pId--)
		{
			// Skip the start of the bifurcation unless at the first bifurcation.
			if(it->first > treeInfo.begin()->first && pId == 0)
			{
				continue;
			}

			vtkMath::Normalize(c0);

			GetDirectionVector(it->first, pId, v1);
			vtkMath::MultiplyScalar(v1, -1);

			DoubleCross(v0, c0, v1, c1);
			vtkMath::Normalize(c1);

			vtkIdType globalPointId = GetGlobalPointId(it->first, pId);

			// Add the twist.
			if(it->first > treeInfo.begin()->first)
			{
				double angleInc = twistAngle / (double)pId;
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->RotateWXYZ(twistDirection * angleInc, v1);
				transform->TransformPoint(c1, c1);
				twistAngle -= angleInc;
			}

			// Store radius vector at point pId.
			radiiVectors->SetTuple(globalPointId, c1);

			// Save v1 as v0 for next iteration.
			vtkMath::Add(v1, zero, v0);

			// Save c1 as c0 for next iteration.
			vtkMath::Add(c1, zero, c0);
		}
	}

#if PRINT_DEBUG
	// Print the map.
	for(std::map<vtkIdType, vtkIdType>::iterator it = lineIdToAverageVectorMap.begin(); it != lineIdToAverageVectorMap.end(); ++it)
	{
		std::cout << it->first << " => " << it->second << std::endl;
		double tuple[3];
		bifurcationVectors->GetTuple(it->second, tuple);
		PrintPoint(tuple); std::cout << std::endl;
	}
#endif

	// Process the terminals.
	for(std::map<vtkIdType, std::vector<vtkIdType> >::iterator it = treeInfo.begin(); it != treeInfo.end(); ++it)
	{
		// Skip all bifurcations.
		if(it->second.size() > 0)
		{
			continue;
		}

		// Process the rest of the line.
		vtkSmartPointer<vtkIdList> lineIds = GetLineIds(it->first);

		vtkIdType parentBifId = GetGlobalPointId(it->first, d_start);

		GetDirectionVector(it->first, d_start, v0);
		vtkMath::MultiplyScalar(v0, -1);

		// c0 from parent bifurcation.
		radiiVectors->GetTuple(parentBifId, c0);

		// Traversing forward.
		for(vtkIdType pId = 1; pId < lineIds->GetNumberOfIds(); pId++)
		{
			vtkMath::Normalize(c0);

			GetDirectionVector(it->first, pId, v1);
			vtkMath::MultiplyScalar(v1, -1);

			vtkIdType globalPointId = GetGlobalPointId(it->first, pId);

			DoubleCross(v0, c0, v1, c1);
			vtkMath::Normalize(c1);

			// Set the length of c1 to the radius value at this point.
			// ==> vtkMath::MultiplyScalar(c1, resampledVesselCentreline->GetPointData()->GetScalars()->GetTuple1(globalPointId));

			// Store radius vector at point pId.
			radiiVectors->SetTuple(globalPointId, c1);

			// Save v1 as v0 for next iteration.
			vtkMath::Add(v1, zero, v0);

			// Save c1 as c0 for next iteration.
			vtkMath::Add(c1, zero, c0);
		}
	}

	// Set radii vectors from normalised vectors to correct magnitude.

	// Root point.
	radiiVectors->GetTuple(0, c0);
	vtkMath::MultiplyScalar(c0, input->GetPointData()->GetScalars()->GetTuple1(0));
	radiiVectors->SetTuple(0, c0);

	lines->InitTraversal();

	vtkSmartPointer<vtkIdList> line = vtkSmartPointer<vtkIdList>::New();

	for(vtkIdType id = 0; id < lines->GetNumberOfCells(); id++)
	{
		lines->GetNextCell(line);

		for(vtkIdType posInLine = 1; posInLine < line->GetNumberOfIds(); posInLine++)
		{
			vtkIdType pId = line->GetId(posInLine);

			radiiVectors->GetTuple(pId, c0);

			// Interpolation of radii vectors.
			if(posInLine < line->GetNumberOfIds() - 1)
			{
				GetDirectionVector(id, posInLine, v0);
				GetDirectionVector(id, posInLine + 1, v1);

				DoubleCross(v0, c0, v1, c1);

				vtkMath::Normalize(c1);
				vtkMath::Add(c0, c1, c0);
				vtkMath::MultiplyScalar(c0, 0.5);

				vtkMath::Normalize(c0);
			}

			vtkMath::MultiplyScalar(c0, input->GetPointData()->GetScalars()->GetTuple1(pId));
			radiiVectors->SetTuple(pId, c0);
		}
	}


	// Shallow copy is good enough?
	output->ShallowCopy(input);

	// TODO: Review whether this is OK to modify/augment input data with vectors. Perhaps the old radii should be stripped from this data and maybe it is cleaner to create new polydata.
	output->GetPointData()->SetVectors(radiiVectors);

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

// Get the vector indicating the direction of the centreline at this location.
// TODO: Perhaps this method should call the other one that overrides it.
void vtkScalarRadiiToVectorsFilter::GetDirectionVector(vtkIdType lineId, LocationType location, double *vector)
{
	vtkSmartPointer<vtkIdList> line = GetLineIds(lineId);

	double p0[3];
	double p1[3];

	if(line->GetNumberOfIds() == 0)
	{
		vtkErrorMacro("Line " << lineId << " is empty.");
	}

	if(location == d_start)
	{
		input->GetPoint(line->GetId(1), p0);
		input->GetPoint(line->GetId(0), p1);
	}
	else if(location == d_end)
	{
		input->GetPoint(line->GetId(line->GetNumberOfIds() - 1), p0);
		input->GetPoint(line->GetId(line->GetNumberOfIds() - 2), p1);
	}
	vtkMath::Subtract(p1, p0, vector);
}

// Get the vector indicating the direction of the centreline at this location.
void vtkScalarRadiiToVectorsFilter::GetDirectionVector(vtkIdType lineId, vtkIdType pointId, double *vector)
{
	vtkSmartPointer<vtkIdList> line = GetLineIds(lineId);

	// TODO: Use error macro to check there are enough points for this.
	assert(pointId >= 0);
	assert(line->GetNumberOfIds() > pointId);

	double p0[3];
	double p1[3];

	// Mapping from a local id in the line to polydata global point id.
	if(pointId != 0)
	{
		input->GetPoint(line->GetId(pointId), p0);
		input->GetPoint(line->GetId(pointId - 1), p1);
		vtkMath::Subtract(p1, p0, vector);
	}
	else
	{
		// This is a bit fishy, because the direction at the fist point is calculated an the direction at the next point.
		GetDirectionVector(lineId, pointId + 1, vector);
	}
}

// Get global point id.
// TODO: Perhaps this method should call the other one that overrides it.
vtkIdType vtkScalarRadiiToVectorsFilter::GetGlobalPointId(vtkIdType lineId, LocationType location)
{
	vtkSmartPointer<vtkIdList> line = GetLineIds(lineId);

	if(line->GetNumberOfIds() == 0)
	{
		vtkErrorMacro("Line " << lineId << " is empty.");
	}

	vtkIdType pId = -1;
	if(location == d_start)
	{
		pId = line->GetId(0);
	}
	else if(location == d_end)
	{
		pId = line->GetId(line->GetNumberOfIds() - 1);
	}
	return pId;
}

// Get global point id.
vtkIdType vtkScalarRadiiToVectorsFilter::GetGlobalPointId(vtkIdType lineId, vtkIdType pointId)
{
	vtkSmartPointer<vtkIdList> line = GetLineIds(lineId);

	if(line->GetNumberOfIds() <= pointId)
	{
		vtkErrorMacro("Point id " << pointId << " is out of bounds.");
	}

	// Mapping from a local id in the line to polydata global point id.
	return line->GetId(pointId);
}

// Get ids for points in this line.
vtkSmartPointer<vtkIdList> vtkScalarRadiiToVectorsFilter::GetLineIds(vtkIdType lineId)
{
	bool lineFound = false;
	vtkSmartPointer<vtkCellArray> lines = input->GetLines();
	vtkSmartPointer<vtkIdList> line = vtkSmartPointer<vtkIdList>::New();
	lines->InitTraversal();
	for(vtkIdType id = 0; id < lines->GetNumberOfCells(); id++)
	{
		lines->GetNextCell(line);
		if(id == lineId)
		{
			lineFound = true;
			break;
		}
	}
	if(lineFound == false)
	{
		vtkErrorMacro("Line " << lineId << " not found.");
	}
	return line;
}
