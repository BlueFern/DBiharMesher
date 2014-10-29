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

	// input = GetInput();

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

	// Obtain tree structure (bifurcations and terminals) from the centreline points, lines.
	for(vtkIdType lineId = 0; lineId < lines->GetNumberOfCells(); lineId++)
	{
		vtkSmartPointer<vtkIdList> lineIds = vtkSmartPointer<vtkIdList>::New(); // TODO: Take this out of the loop.

		lines->GetNextCell(lineIds);
		std::cout << lineId << ": " << lineIds->GetNumberOfIds() << std::endl;

		vtkSmartPointer<vtkIdList> lastId = vtkSmartPointer<vtkIdList>::New();
		lastId->InsertNextId(lineIds->GetId(lineIds->GetNumberOfIds() - 1));

		// WARNING: The code in the GetCellNeighbors method changes the traversal location in the vtkCellArray. In my view it is utterly retarded.
		vtkIdType traverseLocation = input->GetLines()->GetTraversalLocation();

		vtkSmartPointer<vtkIdList> neighbourCellIds = vtkSmartPointer<vtkIdList>::New();
		input->GetCellNeighbors(lineId, lastId, neighbourCellIds);

		// Restore traversal location.
		input->GetLines()->SetTraversalLocation(traverseLocation);

		//std::cout << "Line " << lineId << " shares last point with " << neighbourCellIds->GetNumberOfIds() << " cells." << std::endl;

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

	// TODO: Check this can be removed safely.
	// For now, populate it just in case some tuples get missed.
	for(int i = 0; i < radiiVectors->GetNumberOfTuples(); i++)
	{
		radiiVectors->SetTuple3(i, 0, 0, 0);
	}

	double zero[3] = {0};
	double v0[3];
	double v1[3];
	double c0[3];
	double c1[3];

	// Process the bifurcations to find c0 at each bifurcation point.
	for(std::map<vtkIdType, std::vector<vtkIdType> >::iterator it = treeInfo.begin(); it != treeInfo.end(); ++it)
	{
		// Skip all terminals
		if(it->second.size() == 0)
		{
			continue;
		}

		// Find radius vector c0 at this bifurcation.
		GetDirectionVector(it->second[0], 0, v0);
		vtkMath::MultiplyScalar(v0, -1);
		GetDirectionVector(it->second[1], 0, v1);
		vtkMath::MultiplyScalar(v1, -1);

		// c0 orthogonal to v0 and v1.
		vtkMath::Cross(v0, v1, c0);
		vtkMath::Normalize(c0);

		// Global point id.
		vtkSmartPointer<vtkIdList> lineIds = GetLineIds(it->first);

		// Store bifurcation vector c0.
		radiiVectors->SetTuple(lineIds->GetId(lineIds->GetNumberOfIds() - 1), c0);

		// v0 here is an average of the two vectors.
		vtkMath::Add(v0, v1, v0);
		vtkMath::MultiplyScalar(v0, 0.5);
		vtkMath::Normalize(v0);
		avrgVectors[it->first] = vtkVector3d(v0);
	}

	// First pass to calculate twist angles.
	std::map<vtkIdType, double> twistAngles;
	for(std::map<vtkIdType, std::vector<vtkIdType> >::iterator it = treeInfo.begin(); it != treeInfo.end(); ++it)
	{
		// Skip all terminals
		if(it->second.size() == 0)
		{
			continue;
		}

		vtkSmartPointer<vtkIdList> lineIds = GetLineIds(it->first);

		GetDirectionVector(it->first, lineIds->GetNumberOfIds() - 1, v0);
		vtkMath::MultiplyScalar(v0, -1);

		radiiVectors->GetTuple(lineIds->GetId(lineIds->GetNumberOfIds() - 1), c0);

		// Traversing backwards.
		for(vtkIdType pId = lineIds->GetNumberOfIds() - 2; pId >= 0; pId--)
		{
			// vtkMath::Normalize(c0);

			// Process first point only for inlet.
			if(pId == 0 && it->first != treeInfo.begin()->first)
			{
				break;
			}

			GetDirectionVector(it->first, pId, v1);

			DoubleCross(v0, c0, v1, c1);
			vtkMath::Normalize(c1);

			// Store radius vector at point pId. This is only inlet.
			if(pId == 0 && it->first == treeInfo.begin()->first)
			{
				radiiVectors->SetTuple(lineIds->GetId(pId), c1);
			}

			// Calculate twist angle.
			if(pId == 1 && it->first != treeInfo.begin()->first)
			{
				radiiVectors->GetTuple(lineIds->GetId(0), c0);
				GetDirectionVector(it->first, 0, v0);
				double tmp[3];
				DoubleCross(v0, c0, v1, tmp);

				twistAngles[it->first] = vtkMath::DegreesFromRadians(AngleBetweenVectors(tmp, c1));

				// Figure out twist direction.
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->RotateWXYZ(twistAngles[it->first], v1);
				transform->TransformPoint(tmp, tmp);

				if(vtkMath::DegreesFromRadians(AngleBetweenVectors(tmp, c1)) < angleTolerance)
				{
					twistAngles[it->first] *= -1.0;
				}
				std::cout << twistAngles[it->first] << std::endl;
			}

			// Save v1 as v0 for next iteration.
			vtkMath::Add(v1, zero, v0);

			// Save c1 as c0 for next iteration.
			vtkMath::Add(c1, zero, c0);
		}
	}

	// Second pass to insert radii for all locations other than bifurcations and inlet point.
	for(std::map<vtkIdType, std::vector<vtkIdType> >::iterator it = treeInfo.begin(); it != treeInfo.end(); ++it)
	{
		// Process the rest of the line.
		vtkSmartPointer<vtkIdList> lineIds = GetLineIds(it->first);

		GetDirectionVector(it->first, 0, v0);
		vtkMath::MultiplyScalar(v0, -1);

		// c0 from parent bifurcation.
		radiiVectors->GetTuple(lineIds->GetId(0), c0);

		bool terminal = (it->second.size() == 0);

		// Calculate twist increment.
		double twistAngle = 0.0;
		if(!terminal && (it->first > treeInfo.begin()->first))
		{
			twistAngle = twistAngles[it->first];
		}

		// Traversing forward.
		for(vtkIdType pId = 1; pId < lineIds->GetNumberOfIds(); pId++)
		{
			// If not on terminal, do not insert/overwrite last point.
			if(!terminal && pId == lineIds->GetNumberOfIds() - 1)
			{
				break;
			}

			//vtkMath::Normalize(c0);

			GetDirectionVector(it->first, pId, v1);
			vtkMath::MultiplyScalar(v1, -1);

			DoubleCross(v0, c0, v1, c1);
			vtkMath::Normalize(c1);

			// Add the twist.
			if(!terminal && (it->first > treeInfo.begin()->first))
			{
				double angleInc = twistAngle / (double)(lineIds->GetNumberOfIds() - 1 - pId);
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->RotateWXYZ(angleInc, v1);
				transform->TransformPoint(c1, c1);
				twistAngle -= angleInc;
			}

			// Store radius vector at point pId.
			radiiVectors->SetTuple(lineIds->GetId(pId), c1);

			// Save v1 as v0 for next iteration.
			vtkMath::Add(v1, zero, v0);

			// Save c1 as c0 for next iteration.
			vtkMath::Add(c1, zero, c0);
		}
	}

	// Root point.
	radiiVectors->GetTuple(0, c0);
	vtkMath::MultiplyScalar(c0, input->GetPointData()->GetScalars()->GetTuple1(0));
	radiiVectors->SetTuple(0, c0);

	// Third pass to interpolate radii on sharp angles and set correct magnitude.
	for(std::map<vtkIdType, std::vector<vtkIdType> >::iterator it = treeInfo.begin(); it != treeInfo.end(); ++it)
	{
		// Process the rest of the line.
		vtkSmartPointer<vtkIdList> lineIds = GetLineIds(it->first);

		for(vtkIdType pId = 1; pId < lineIds->GetNumberOfIds(); pId++)
		{
			radiiVectors->GetTuple(lineIds->GetId(pId), c0);

			// Interpolation of radii vectors.
			if(pId < lineIds->GetNumberOfIds() - 2)
			{
				GetDirectionVector(it->first, pId, v0);
				GetDirectionVector(it->first, pId + 1, v1);

				DoubleCross(v0, c0, v1, c1);

				vtkMath::Normalize(c1);
				vtkMath::Add(c0, c1, c0);
				vtkMath::MultiplyScalar(c0, 0.5);

				vtkMath::Normalize(c0);
			}

			vtkMath::MultiplyScalar(c0, input->GetPointData()->GetScalars()->GetTuple1(lineIds->GetId(pId)));
			radiiVectors->SetTuple(lineIds->GetId(pId), c0);
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
void vtkScalarRadiiToVectorsFilter::GetDirectionVector(vtkIdType lineId, vtkIdType pointId, double *vector)
{
	vtkSmartPointer<vtkIdList> line = GetLineIds(lineId);

	if(pointId < 0 || pointId > line->GetNumberOfIds())
	{
		vtkErrorMacro("Point id " << pointId << " is out of bounds for line id " << lineId << ".");
	}

	if(pointId == line->GetNumberOfIds() - 1 && avrgVectors.find(lineId) != avrgVectors.end())
	{
		vector[0] = avrgVectors[lineId][0];
		vector[1] = avrgVectors[lineId][1];
		vector[2] = avrgVectors[lineId][2];
		return;
	}

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
