/**
 * Program: vtkEndPointIdsToDbiharPatchFilter.
 */

#include <vector>
#include <algorithm>

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkCallbackCommand.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>

#include <vtkAppendPolyData.h>
#include <vtkDataObject.h>
#include <vtkPointData.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkTransform.h>
#include <vtkIdList.h>
#include <vtkMath.h>

#include "vtkEndPointIdsToDbiharPatchFilter.h"
#include "showPolyData.h"

#define PRINT_DEBUG 0

// TODO: Move to a lib. Code is duplicated.
void DoubleCross1(const double v0[3], const double c0[3], const double v1[3], double c1[3])
{
	vtkMath::Cross(c0, v0, c1);
	vtkMath::Cross(v1, c1, c1);
}

vtkStandardNewMacro(vtkEndPointIdsToDbiharPatchFilter);

// TODO: Declare this in a better place shared between classes.
const char *vtkEndPointIdsToDbiharPatchFilter::RADII_ARR_NAME = {"radiiVectors"};

vtkEndPointIdsToDbiharPatchFilter::vtkEndPointIdsToDbiharPatchFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	// this->SegmentIdList = vtkSmartPointer<vtkIdList>::New();

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkEndPointIdsToDbiharPatchFilter::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	if(EndPointIdsList->GetNumberOfIds() == 0)
	{
		vtkErrorMacro("Segment id list is empty.");
	}

	if((NumberOfRadialQuads & 1) != 0)
	{
		vtkErrorMacro("Number of radial quads must be even.")
	}

	if(NumberOfRadialQuads < 4)
	{
		vtkErrorMacro("Number of radial quads must be at least 4.")
	}

	// Get the input and output.
	input = vtkPolyData::GetData(inputVector[0], 0);
	input->BuildLinks();

	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

#if 0
	vtkSmartPointer<vtkCellArray> lines = input->GetLines();
	lines->InitTraversal();

	// Obtain tree structure (bifurcations and terminals) from the centreline points, lines.
	for(vtkIdType lineId = 0; lineId < lines->GetNumberOfCells(); lineId++)
	{
		vtkSmartPointer<vtkIdList> lineIds = vtkSmartPointer<vtkIdList>::New(); // TODO: Take this out of the loop.

		lines->GetNextCell(lineIds);
		// std::cout << lineId << ": " << lineIds->GetNumberOfIds() << std::endl;

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
#endif

	if(EndPointIdsList->GetNumberOfIds() < 2)
	{
		vtkErrorMacro("End point id list is expected to have at least two points.");
	}
	else if(EndPointIdsList->GetNumberOfIds() > 3)
	{
		vtkWarningMacro("Patch generation for more than two branches at bifurcations has not been tested.");
	}

	std::vector<std::vector<vtkIdType> > branchIds;

	// Process segment id list to produce patches.
	for(int i = 0; i < EndPointIdsList->GetNumberOfIds(); i++)
	{
		branchIds.push_back(std::vector<vtkIdType>());
	}

	// TODO: Need to differentiate between straight segments and bifurcations.

	// 1. Find current cell id.
	// Temporary storage for cell ids shared by current point.
	vtkSmartPointer<vtkIdList> tmpIdList = vtkSmartPointer<vtkIdList>::New();
	input->GetPointCells(EndPointIdsList->GetId(0), tmpIdList);
	if(tmpIdList->GetNumberOfIds() > 1)
	{
		vtkErrorMacro("First point id is shared between multiple cells. Dbihar patch must not start on a bifurcation.");
	}
	vtkIdType trunkCellId = tmpIdList->GetId(0);

	// 2. Find last point id for the nearest downstream bifurcation.
	tmpIdList = input->GetCell(trunkCellId)->GetPointIds();
	vtkIdType bifurcationId = tmpIdList->GetId(tmpIdList->GetNumberOfIds() - 1);

	// 3. For trunkCellId traverse back to collect ids. Exclude bifurcation id.
	for(int id = tmpIdList->GetNumberOfIds() - 2; id >= EndPointIdsList->GetId(0); id--)
	{
		branchIds[0].push_back(tmpIdList->GetId(id));
	}
	std::reverse(branchIds[0].begin(), branchIds[0].end());

	// 4. Find branch ids.
	vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
	ptIds->InsertNextId(bifurcationId);
	input->GetCellNeighbors(trunkCellId, ptIds, cellIds);

	// 5. Collect branch ids.
	for(int spineId = 1; spineId < branchIds.size(); spineId++)
	{
		tmpIdList = input->GetCell(cellIds->GetId(spineId - 1))->GetPointIds();
		vtkIdType endId = EndPointIdsList->GetId(spineId);
		for(int id = 1; id < tmpIdList->GetNumberOfIds(); id++)
		{
			vtkIdType pId = tmpIdList->GetId(id);
			branchIds[spineId].push_back(pId);
			if(pId == endId)
			{
				break;
			}
		}
	}

#if 0
	for(int spineId = 0; spineId < branchIds.size(); spineId++)
	{
		for(std::vector<vtkIdType>::iterator it = branchIds[spineId].begin(); it != branchIds[spineId].end(); ++it)
		{
			std::cout << *it << " ";
		}
		std::cout << std::endl;
	}
#endif

	// Now assemble branch ids.
	for(int i = 0; i < EndPointIdsList->GetNumberOfIds(); i++)
	{
		spineIds.push_back(std::vector<vtkIdType>());
	}

	// First patch id list.
	spineIds[0].insert(spineIds[0].end(), branchIds[0].begin(), branchIds[0].end());
	spineIds[0].push_back(bifurcationId);
	spineIds[0].insert(spineIds[0].end(), branchIds[1].begin(), branchIds[1].end());

	// Last patch id list.
	spineIds[spineIds.size() - 1].insert(spineIds[spineIds.size() - 1].end(), branchIds[branchIds.size() - 1].rbegin(), branchIds[branchIds.size() - 1].rend());
	spineIds[spineIds.size() - 1].push_back(bifurcationId);
	spineIds[spineIds.size() - 1].insert(spineIds[spineIds.size() - 1].end(), branchIds[0].rbegin(), branchIds[0].rend());

	// Bifurcation saddle patch id list.
	// Looping is not required here, but with the use of a loop the code can be extended to more than two branches.
	for(int id = 1; id < cellIds->GetNumberOfIds(); id++)
	{
		spineIds[id].insert(spineIds[id].end(), branchIds[id].rbegin(), branchIds[id].rend());
		spineIds[id].push_back(bifurcationId);
		spineIds[id].insert(spineIds[id].end(), branchIds[id + 1].begin(), branchIds[id + 1].end());
	}

#if 1
	for(int spineId = 0; spineId < spineIds.size(); spineId++)
	{
		for(std::vector<vtkIdType>::iterator it = spineIds[spineId].begin(); it != spineIds[spineId].end(); ++it)
		{
			std::cout << *it << " ";
		}
		std::cout << std::endl;
	}
#endif

#if 0
	std::vector<vtkSmartPointer<vtkPolyData> > inputPatches;
	for(int spineId = 0; spineId < spineIds.size(); spineId++)
	{
		inputPatches[spineId].Take(vtkPolyData::New());
	}
#endif

	vtkSmartPointer<vtkAppendPolyData> appendPolyDataFilter = vtkSmartPointer<vtkAppendPolyData>::New();

	for(int spineId = 0; spineId < spineIds.size(); spineId++)
	{
		// TODO: It seems to make sense to have this code for patch construction in a stand-alone filter.

		vtkSmartPointer<vtkPoints> patchPoints = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPolyLine> patchBoundary = vtkSmartPointer<vtkPolyLine>::New();
		vtkSmartPointer<vtkIdList> verts = vtkSmartPointer<vtkIdList>::New();

		int spineLength = spineIds[spineId].size();

		int numPtIds = spineLength * 2 + this->NumberOfRadialQuads * 2 - 2;

		for(vtkIdType ptId = 0, spinePtId = 0; ptId < numPtIds; ptId++)
		{
			// Point is declared inside the loop to make sure it is (0,0,0) at the start of every iteration.
			double point[3] = {0.0};
			// Derivative is declared inside the loop to make sure it is (0,0,0) at the start of every iteration.
			double deriv[3] = {0.0};

			if(ptId < this->NumberOfRadialQuads) // Number of points along this edge is number of quads + 1.
			{
				vtkIdType localId = ptId;
				std::cout << "LC: " << localId << std::endl;

				// Inserting the first arc.
				double parametricCoord = ptId / (double)this->NumberOfRadialQuads;

				double p0[3];
				double p1[3];
				double v0[3];
				double c0[3];

				// Translation comes from the first point.
				input->GetPoint(spineIds[spineId][0], p0);

				// Axis of rotation comes from the centreline direction at the first point.
				input->GetPoint(spineIds[spineId][1], p1);
				vtkMath::Subtract(p1, p0, v0);

				// Get the first radius.
				double radius[3];
				vtkSmartPointer<vtkDoubleArray> radiiArray = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetVectors(RADII_ARR_NAME));
				radiiArray->GetTuple(spineIds[spineId][0], radius);

				// Get the rotation axis.
				DoubleCross1(radius, v0, radius, c0);

				// Angle of rotation comes from the parametric coordinate along the arc.
				double angle = vtkMath::Pi() * parametricCoord;

				// Assemble local transform.
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->Translate(p0);
				transform->RotateWXYZ(vtkMath::DegreesFromRadians(angle), c0);

				// Transform the radius.
				transform->TransformPoint(radius, point);
			}
			else if(ptId < this->NumberOfRadialQuads + spineLength - 1)
			{
				// Inserting left seam.
				vtkIdType localId = ptId - this->NumberOfRadialQuads;
				std::cout << "LS: " << localId << std::endl;

				double p0[3];

				// Translation comes from the first point.
				input->GetPoint(spineIds[spineId][localId], p0);
				std::cout << "--> " << spineIds[spineId][localId] << std::endl;

				// Get the first radius.
				double radius[3];
				vtkSmartPointer<vtkDoubleArray> radiiArray = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetVectors(RADII_ARR_NAME));
				radiiArray->GetTuple(spineIds[spineId][localId], radius);

				// Flip the radius.
				vtkMath::MultiplyScalar(radius, -1.0);

				// Translate the radius.
				vtkMath::Add(p0, radius, point);
			}
			else if(ptId < this->NumberOfRadialQuads * 2 + spineLength - 1)
			{
				// Inserting the second arc.
				vtkIdType localId = ptId - this->NumberOfRadialQuads - (spineLength - 1);
				std::cout << "UC: " << localId << std::endl;

				double parametricCoord = localId / (double)this->NumberOfRadialQuads;
				// Inserting the first arc.

				double p0[3];
				double p1[3];
				double v0[3];
				double c0[3];

				// Translation comes from the first point.
				input->GetPoint(spineIds[spineId][spineLength - 1], p0);

				std::cout << "--> " << spineIds[spineId][spineLength - 1] << std::endl;

				// Axis of rotation comes from the centreline direction at the first point.
				input->GetPoint(spineIds[spineId][spineLength - 2], p1);
				vtkMath::Subtract(p1, p0, v0);

				// Get the last radius.
				double radius[3];
				vtkSmartPointer<vtkDoubleArray> radiiArray = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetVectors(RADII_ARR_NAME));
				radiiArray->GetTuple(spineIds[spineId][spineLength - 1], radius);
				vtkMath::MultiplyScalar(radius, -1.0);

				// Get the rotation axis.
				DoubleCross1(radius, v0, radius, c0);

				// Angle of rotation comes from the parametric coordinate along the arc.
				double angle = vtkMath::Pi() * parametricCoord;

				// Assemble local transform.
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->Translate(p0);
				transform->RotateWXYZ(vtkMath::DegreesFromRadians(angle), c0);

				// Transform the radius.
				transform->TransformPoint(radius, point);
			}
			else
			{
				// Inserting the last seam.
				vtkIdType localId = ptId - this->NumberOfRadialQuads * 2  - (spineLength - 1);
				localId = std::fabs(localId - (spineLength - 1));
				std::cout << "RS: " << localId << std::endl;

				double p0[3];

				// Translation comes from the first point.
				input->GetPoint(spineIds[spineId][localId], p0);
				std::cout << "--> " << spineIds[spineId][localId] << std::endl;

				// Get the first radius.
				double radius[3];
				vtkSmartPointer<vtkDoubleArray> radiiArray = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetVectors(RADII_ARR_NAME));
				radiiArray->GetTuple(spineIds[spineId][localId], radius);

				// Translate the radius.
				vtkMath::Add(p0, radius, point);
			}

			vtkIdType id = patchPoints->InsertNextPoint(point);
			// Sanity check.
			assert(id == ptId);
			patchBoundary->GetPointIds()->InsertNextId(ptId);
			verts->InsertNextId(ptId);
		}
		patchBoundary->GetPointIds()->InsertNextId(0);

		vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
		boundaries->InsertNextCell(patchBoundary);

		vtkSmartPointer<vtkCellArray> vertsArray = vtkSmartPointer<vtkCellArray>::New();
		vertsArray->InsertNextCell(verts);

		vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
		inputPatch->SetPoints(patchPoints);
		inputPatch->SetLines(boundaries);
		inputPatch->SetVerts(vertsArray);

		std::cout << inputPatch->GetNumberOfPoints() << " =?= " << numPtIds << std::endl;

		showPolyData(inputPatch, NULL);

		appendPolyDataFilter->AddInputData(inputPatch);

		// TODO: Insert the derivatives.

		//inputPatches[spineId]->DeepCopy(inputPatch);
	}

	appendPolyDataFilter->Update();
	showPolyData(appendPolyDataFilter->GetOutput(), NULL);

	// TODO: Use vtkDbiharPatchFiler to obtain output patches.

	// TODO: Merge output patches into one vtkPolyData object.

	// TODO: Implement progress updates.
	// this->UpdateProgress(static_cast<double>(pId)/static_cast<double>(numPIds));

	output->DeepCopy(output);

	return 1;
}

void vtkEndPointIdsToDbiharPatchFilter::PrintSelf(ostream& os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);

	// TODO: Print idList.
}

void vtkEndPointIdsToDbiharPatchFilter::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkEndPointIdsToDbiharPatchFilter* filter = static_cast<vtkEndPointIdsToDbiharPatchFilter *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}
