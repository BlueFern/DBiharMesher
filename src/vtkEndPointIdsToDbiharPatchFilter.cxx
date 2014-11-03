/**
 * Program: vtkEndPointIdsToDbiharPatchFilter.
 */

#include <vector>
#include <algorithm>

#include <vtkObjectFactory.h>
#include <vtkCallbackCommand.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>

#include <vtkCellArray.h>
#include <vtkIdList.h>

#include "vtkEndPointIdsToDbiharPatchFilter.h"

#define PRINT_DEBUG 0

vtkStandardNewMacro(vtkEndPointIdsToDbiharPatchFilter);

#if 0
const char *vtkEndPointIdsToDbiharPatchFilter::RADII_ARR_NAME = {"radiiVectors"};
#endif

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

	if(EndPointIdsList->GetNumberOfIds() > 2)
	{
		vtkErrorMacro("End point id list is expected to have between two and three points");
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
	spineIds[spineIds.size() - 1].insert(spineIds[spineIds.size() - 1].end(), branchIds[0].begin(), branchIds[0].end());
	spineIds[spineIds.size() - 1].push_back(bifurcationId);
	spineIds[spineIds.size() - 1].insert(spineIds[spineIds.size() - 1].end(), branchIds[branchIds.size() - 1].begin(), branchIds[branchIds.size() - 1].end());

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

	std::vector<vtkSmartPointer<vtkPolyData> > inputPatches;
	for(int spineId = 0; spineId < spineIds.size(); spineId++)
	{
		inputPatches[spineId] = vtkSmartPointer<vtkPolyData>::New();
	}

	// TODO: Construct input patches.

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
