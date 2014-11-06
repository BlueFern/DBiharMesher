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
#include <vtkStructuredGrid.h>
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
#include "vtkDbiharPatchFilter.h"
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
const char *vtkEndPointIdsToDbiharPatchFilter::DERIV_ARR_NAME = {"derivVectors"};

vtkEndPointIdsToDbiharPatchFilter::vtkEndPointIdsToDbiharPatchFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	// this->SegmentIdList = vtkSmartPointer<vtkIdList>::New();
	this->cEdgeScaling = 3;
	this->yEdgeScaling = 6.0;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkEndPointIdsToDbiharPatchFilter::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	bool bifurcation = true;
	if(EndPointIdsList->GetNumberOfIds() == 0)
	{
		vtkErrorMacro("Segment id list is empty.");
	}
	else if(EndPointIdsList->GetNumberOfIds() == 1)
	{
		vtkErrorMacro("Segment id list must have at least two point ids.")
	}
	else if(EndPointIdsList->GetNumberOfIds() == 2)
	{
		bifurcation = false;
	}
	else if(EndPointIdsList->GetNumberOfIds() > 3)
	{
		vtkWarningMacro("Patch generation for more than two branches at bifurcations has not been tested.");
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

	// Process segment id list to produce patches.

	vtkSmartPointer<vtkIdList> tmpIdList = vtkSmartPointer<vtkIdList>::New();

	// 1. Find current cell id.
	input->GetPointCells(EndPointIdsList->GetId(0), tmpIdList);
	if(tmpIdList->GetNumberOfIds() > 1)
	{
		vtkErrorMacro("First point id is shared between multiple cells. Dbihar patch must not start on a bifurcation.");
	}
	vtkIdType trunkCellId = tmpIdList->GetId(0);

	// All ids in the trunk cell.
	vtkSmartPointer<vtkIdList> trunkIdCellIdList = input->GetCell(trunkCellId)->GetPointIds();

	// Now assemble branch ids.
	for(int i = 0; i < EndPointIdsList->GetNumberOfIds(); i++)
	{
		spineIds.push_back(std::vector<vtkIdType>());
	}

	if(!bifurcation)
	{
		// Collecting ids for a straight segment.

		vtkIdType locId = -1;
		vtkIdType ptId = -1;

		// Find the location of the first point.
		while(true)
		{
			// TODO: Verify we are not going to fall off the end of the list and not land in a bifurcation.
			ptId = trunkIdCellIdList->GetId(++locId);
			if(ptId == EndPointIdsList->GetId(0))
			{
				break;
			}
		}

		// Remember the point ids for this trunk.
		while(true)
		{
			// TODO: Verify we are not going to fall off the end of the list and not land in a bifurcation.
			spineIds[0].push_back(ptId);
			spineIds[1].push_back(ptId);
			if(ptId == EndPointIdsList->GetId(1))
			{
				break;
			}
			else
			{
				ptId = trunkIdCellIdList->GetId(++locId);
			}
		}
		// For the second patch the ids need to be reversed.
		std::reverse(spineIds[1].rbegin(), spineIds[1].rend());
	}
	else
	{
		// Collecting ids for a bifurcation.

		std::vector<std::vector<vtkIdType> > branchIds;
		for(int i = 0; i < EndPointIdsList->GetNumberOfIds(); i++)
		{
			branchIds.push_back(std::vector<vtkIdType>());
		}

		// TODO: This code is not robust, because we are assuming the end of the list is the one bifurcation. Find a better way!
		// For trunkCellId traverse back to collect ids. Exclude bifurcation id.
		for(int id = trunkIdCellIdList->GetNumberOfIds() - 2; id >= EndPointIdsList->GetId(0); id--)
		{
			branchIds[0].push_back(trunkIdCellIdList->GetId(id));
		}
		std::reverse(branchIds[0].begin(), branchIds[0].end());

		// TODO: This code is not robust, because we are assuming the end of the list is the one bifurcation. Find a better way!
		// Find last point id for the nearest downstream bifurcation.
		vtkIdType bifurcationId = trunkIdCellIdList->GetId(trunkIdCellIdList->GetNumberOfIds() - 1);

		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
		ptIds->InsertNextId(bifurcationId);
		input->GetCellNeighbors(trunkCellId, ptIds, cellIds);

		// Collect branch ids.
		for(int spineId = 1; spineId < branchIds.size(); spineId++)
		{
			// TODO: Need error checking to catch us falling of the end of the id list.
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

		// First patch id list.
		spineIds[0].insert(spineIds[0].end(), branchIds[0].begin(), branchIds[0].end());
		spineIds[0].push_back(bifurcationId);
		spineIds[0].insert(spineIds[0].end(), branchIds[1].begin(), branchIds[1].end());

		// Last patch id list.
		spineIds[spineIds.size() - 1].insert(spineIds[spineIds.size() - 1].end(), branchIds[branchIds.size() - 1].rbegin(), branchIds[branchIds.size() - 1].rend());
		spineIds[spineIds.size() - 1].push_back(bifurcationId);
		spineIds[spineIds.size() - 1].insert(spineIds[spineIds.size() - 1].end(), branchIds[0].rbegin(), branchIds[0].rend());

		// Bifurcation saddle patch(es) id list.
		// Looping is not required here, but with the use of a loop the code can be extended to more than two branches.
		for(int id = 1; id < cellIds->GetNumberOfIds(); id++)
		{
			spineIds[id].insert(spineIds[id].end(), branchIds[id].rbegin(), branchIds[id].rend());
			spineIds[id].push_back(bifurcationId);
			spineIds[id].insert(spineIds[id].end(), branchIds[id + 1].begin(), branchIds[id + 1].end());
		}
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
	//vtkSmartPointer<vtkStructuredGridAppend> appendStructuredGridFilter = vtkSmartPointer<vtkStructuredGridAppend>::New();

	for(int spineId = 0; spineId < spineIds.size(); spineId++)
	{
		// Points and lines.
		vtkSmartPointer<vtkPoints> patchPoints = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPolyLine> patchBoundary = vtkSmartPointer<vtkPolyLine>::New();
		//vtkSmartPointer<vtkIdList> verts = vtkSmartPointer<vtkIdList>::New();

		// Derivatives.
		vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
		derivatives->SetName(DERIV_ARR_NAME);
		derivatives->SetNumberOfComponents(3);

		int spineLength = spineIds[spineId].size();

		int numPtIds = spineLength * 2 + this->NumberOfRadialQuads * 2 - 2;

		const double zero[3] = {0};
		for(vtkIdType ptId = 0, spinePtId = 0; ptId < numPtIds; ptId++)
		{
			// Point is declared inside the loop to make sure it is (0,0,0) at the start of every iteration.
			double point[3] = {0.0};
			// Derivative is declared inside the loop to make sure it is (0,0,0) at the start of every iteration.
			double deriv[3] = {0.0};

			if(ptId < this->NumberOfRadialQuads) // Number of points along this edge is number of quads + 1.
			{
				// Inserting the first arc.

				vtkIdType localId = ptId;
				std::cout << "LC: " << localId << std::endl;

				double parametricCoord = localId / (double)this->NumberOfRadialQuads;

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

				if(localId != 0)
				{
					// If not at the patch corner point, we need to insert a derivative.

					// It is parallel to the end point direction vector.
					vtkMath::Add(zero, v0, deriv);
					vtkMath::MultiplyScalar(deriv, -1.0);
					vtkMath::Normalize(deriv);
					// Scale derivative vector magnitude.
					// TODO: Vector magnitude should be proportional to the length of the patch?
					vtkMath::MultiplyScalar(deriv, cEdgeScaling);
				}
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

				if(localId != 0)
				{
					// If not at the patch corner point, we need to insert a derivative.
					double p0[3];
					double p1[3];
					double v0[3];

					input->GetPoint(spineIds[spineId][localId], p0);
					input->GetPoint(spineIds[spineId][localId + 1], p1);
					// Get the centreline direction at the current point.
					vtkMath::Subtract(p1, p0, v0);

					// It is perpendicular to the radius vector and the local direction.
					vtkMath::Cross(v0, radius, deriv);
					vtkMath::Normalize(deriv);

					double scaling = vtkMath::Norm(radius) * yEdgeScaling;
					// Scale derivative vector magnitude.
					vtkMath::MultiplyScalar(deriv, scaling);
				}
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

				if(localId != 0)
				{
					// If not at the patch corner point, we need to insert a derivative.
					// It is parallel to the end point direction vector.
					vtkMath::Add(zero, v0, deriv);
					vtkMath::MultiplyScalar(deriv, -1.0);
					vtkMath::Normalize(deriv);
					// Scale derivative vector magnitude.
					// TODO: Vector magnitude should be proportional to the length of the patch?
					vtkMath::MultiplyScalar(deriv, cEdgeScaling);
				}
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

				if(localId != spineLength - 1)
				{
					// If not at the patch corner point, we need to insert a derivative.
					double p0[3];
					double p1[3];
					double v0[3];

					input->GetPoint(spineIds[spineId][localId], p0);
					input->GetPoint(spineIds[spineId][localId + 1], p1);

					// Get the centreline direction at the current point.
					vtkMath::Subtract(p1, p0, v0);

					// It is perpendicular to the radius vector and the local direction.
					vtkMath::Cross(v0, radius, deriv);
					vtkMath::MultiplyScalar(deriv, -1.0);
					vtkMath::Normalize(deriv);

					double scaling = vtkMath::Norm(radius) * yEdgeScaling;
					// Scale derivative vector magnitude.
					vtkMath::MultiplyScalar(deriv, scaling);
				}
			}

			vtkIdType id = patchPoints->InsertNextPoint(point);
			// Sanity check.
			assert(id == ptId);
			patchBoundary->GetPointIds()->InsertNextId(ptId);
			//verts->InsertNextId(ptId);
			derivatives->InsertNextTuple(deriv);
		}
		patchBoundary->GetPointIds()->InsertNextId(0);

		vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
		boundaries->InsertNextCell(patchBoundary);

		//vtkSmartPointer<vtkCellArray> vertsArray = vtkSmartPointer<vtkCellArray>::New();
		//vertsArray->InsertNextCell(verts);

		vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
		inputPatch->SetPoints(patchPoints);
		inputPatch->SetLines(boundaries);
		//inputPatch->SetVerts(vertsArray);
		inputPatch->GetPointData()->SetVectors(derivatives);

		std::cout << inputPatch->GetNumberOfPoints() << " =?= " << numPtIds << std::endl;

		showPolyData(inputPatch, NULL, 1.0);

		vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

		// Set the bounds of the UV space.
		patchFilter->SetA(0.0);
		patchFilter->SetB(2.0/3.0);
		patchFilter->SetC(0.0);
		patchFilter->SetD(vtkMath::Pi());
		// Set the boundary conditions.
		patchFilter->SetMQuads(this->NumberOfRadialQuads);
		patchFilter->SetNQuads(spineLength - 1);
		// Set solution method.
		patchFilter->SetIFlag(2);

		patchFilter->SetInputData(inputPatch);
		patchFilter->Update();

		vtkPolyData *outputPatch = patchFilter->GetOutput();

		vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
		structuredGrid->SetDimensions(this->NumberOfRadialQuads + 1, spineLength, 1);
		structuredGrid->SetPoints(outputPatch->GetPoints());

		showPolyData(inputPatch, structuredGrid);

		appendPolyDataFilter->AddInputData(inputPatch);

		//inputPatches[spineId]->DeepCopy(inputPatch);
	}

	appendPolyDataFilter->Update();
	showPolyData(appendPolyDataFilter->GetOutput(), NULL, 1.0);

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
