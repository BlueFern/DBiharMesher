/**
 * Program: vtkEndPointIdsToDbiharPatchFilter.
 */

#include <vector>
#include <string>
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

//#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
//        ( std::ostringstream() << std::dec << x ) ).str()

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
		exit(EXIT_FAILURE);
	}
	else if(EndPointIdsList->GetNumberOfIds() == 1)
	{
		vtkErrorMacro("Segment id list must have at least two point ids.");
		exit(EXIT_FAILURE);
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
		vtkErrorMacro("Number of radial quads must be even.");
		exit(EXIT_FAILURE);
	}

	if(NumberOfRadialQuads < 4)
	{
		vtkErrorMacro("Number of radial quads must be at least 4.");
		exit(EXIT_FAILURE);
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
		exit(EXIT_FAILURE);
	}
	vtkIdType trunkCellId = tmpIdList->GetId(0);

	// All ids in the trunk cell.
	vtkSmartPointer<vtkIdList> trunkIdCellIdList = input->GetCell(trunkCellId)->GetPointIds();

	// Now assemble branch ids.
	for(int i = 0; i < EndPointIdsList->GetNumberOfIds(); i++)
	{
		spineIds.push_back(std::vector<vtkIdType>());
	}

	vtkIdType bifurcationId = -1;

	if(!bifurcation)
	{
		// Collecting ids for a straight segment.
		vtkIdType locId = -1;
		vtkIdType ptId = -1;
		vtkIdType maxLocId = trunkIdCellIdList->GetNumberOfIds() - 1;

		// Find the location of the first point.
		vtkIdType firstId = EndPointIdsList->GetId(0);
		while(true)
		{
			// TODO: Verify we are not going to fall off the end of the list and not land in a bifurcation.
			++locId;
			if(locId > maxLocId)
			{
				vtkErrorMacro("Unable to find id " << firstId << " in the id list.");
				exit(EXIT_FAILURE);
			}
			ptId = trunkIdCellIdList->GetId(locId);
			if(ptId == firstId)
			{
				break;
			}
		}

		// Remember the point ids for this trunk.
		vtkIdType lastId = EndPointIdsList->GetId(1);
		while(true)
		{
			// TODO: Verify we are not going to fall off the end of the list and not land in a bifurcation.
			spineIds[0].push_back(ptId);
			spineIds[1].push_back(ptId);
			if(ptId == lastId)
			{
				break;
			}
			else
			{
				++locId;
				if(locId > maxLocId)
				{
					vtkErrorMacro("Unable to find id " << lastId << " in the id list.");
					exit(EXIT_FAILURE);
				}
				ptId = trunkIdCellIdList->GetId(locId);
			}
		}

		std::cout << "Reversing list..." << std::endl;

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

		// TODO: This code is not robust, because we are assuming the end of the list is the bifurcation we are interested in. Find a better way!
		// Find last point id for the nearest downstream bifurcation.
		bifurcationId = trunkIdCellIdList->GetId(trunkIdCellIdList->GetNumberOfIds() - 1);

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
	for(int spineId = 0; spineId < spineIds.size(); spineId++)
	{
		inputPatches[spineId].Take(vtkPolyData::New());
	}
#endif

	std::vector<vtkSmartPointer<vtkPolyData> > inputPatches;
	std::vector<vtkSmartPointer<vtkStructuredGrid> > outputGrids;

	vtkSmartPointer<vtkAppendPolyData> appendPolyDataFilter = vtkSmartPointer<vtkAppendPolyData>::New();

	vtkSmartPointer<vtkDoubleArray> radiiArray = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetVectors(RADII_ARR_NAME));

	for(int spineId = 0; spineId < spineIds.size(); spineId++)
	{
		// Points and lines.
		vtkSmartPointer<vtkPoints> patchPoints = vtkSmartPointer<vtkPoints>::New();
		vtkSmartPointer<vtkPolyLine> patchBoundary = vtkSmartPointer<vtkPolyLine>::New();

		// Derivatives.
		vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
		derivatives->SetName(DERIV_ARR_NAME);
		derivatives->SetNumberOfComponents(3);

		int spineLength = spineIds[spineId].size();

		int numPtIds = spineLength * 2 + this->NumberOfRadialQuads * 2 - 2;

		vtkIdType bifurcationPos = -1;
		vtkIdType rightBifurcationDerivId = -1;
		vtkIdType leftBifurcationDerivId = -1;

		// Find position of the bifurcation id in the spineId.
		if(bifurcation)
		{
			for(int ptId = 0; ptId < spineLength; ptId++)
			{
				if(bifurcationId == spineIds[spineId][ptId])
				{
					bifurcationPos = ptId;
					std::cout << "*** " << bifurcationPos << " ***" << std::endl;
					break;
				}
			}
			rightBifurcationDerivId = this->NumberOfRadialQuads + bifurcationPos;
			leftBifurcationDerivId = numPtIds - bifurcationPos;
		}

		const double zero[3] = {0};
		double point[3] = {0.0};

		// Temporary buffers.
		double r[3];
		double p0[3];
		double p1[3];
		double v0[3];
		double c0[3];

		// TODO: Interpolate the derivatives (similar to what we do with radii) between adjacent points. Alternatively, we should try
		// calculating the derivatives only for the centreline and translating them to both edges. They would still have to be interpolated.

		for(vtkIdType ptId = 0, spinePtId = 0; ptId < numPtIds; ptId++)
		{
			// Derivative is declared inside the loop to make sure it is (0,0,0) at the start of every iteration.
			double deriv[3] = {0.0};

			if(ptId < this->NumberOfRadialQuads)
			{
				// Inserting the first (lower edge) arc.

				vtkIdType localId = ptId;
				std::cout << "LC: " << localId << std::endl;

				double parametricCoord = localId / (double)this->NumberOfRadialQuads;

				// Translation comes from the first point.
				input->GetPoint(spineIds[spineId][0], p0);

				// Axis of rotation comes from the centreline direction at the first point.
				input->GetPoint(spineIds[spineId][1], p1);
				vtkMath::Subtract(p1, p0, v0);

				// Get the first radius.
				radiiArray->GetTuple(spineIds[spineId][0], r);

				// Get the rotation axis.
				DoubleCross1(r, v0, r, c0);

				// Angle of rotation comes from the parametric coordinate along the arc.
				double angle = vtkMath::Pi() * parametricCoord;

				// Assemble local transform.
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->Translate(p0);
				transform->RotateWXYZ(vtkMath::DegreesFromRadians(angle), c0);

				// Transform the radius.
				transform->TransformPoint(r, point);

				// If not at the patch corner point, we need to insert a derivative.
				if(localId != 0)
				{
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
				// Inserting along the "left" patch edge.

				vtkIdType localId = ptId - this->NumberOfRadialQuads;
				std::cout << "LS: " << localId << std::endl;

				// Translation comes from the first point.
				input->GetPoint(spineIds[spineId][localId], p0);

				// Get the current radius.
				radiiArray->GetTuple(spineIds[spineId][localId], r);

				// Flip the radius.
				vtkMath::MultiplyScalar(r, -1.0);

				// Translate the radius.
				vtkMath::Add(p0, r, point);

				// If not at the patch corner point, we need to insert a derivative.
				if(localId != 0)
				{
					// Get the centreline direction at the current point.
					input->GetPoint(spineIds[spineId][localId], p0);
					input->GetPoint(spineIds[spineId][localId + 1], p1);
					vtkMath::Subtract(p1, p0, v0);

					// The derivative is perpendicular to the radius vector and the local direction.
					vtkMath::Cross(v0, r, deriv);
					vtkMath::Normalize(deriv);

					double scaling = vtkMath::Norm(r) * yEdgeScaling;
					// Scale derivative vector magnitude.
					vtkMath::MultiplyScalar(deriv, scaling);
				}
			}
			else if(ptId < this->NumberOfRadialQuads * 2 + spineLength - 1)
			{
				// Inserting the second (upper edge) arc.

				vtkIdType localId = ptId - this->NumberOfRadialQuads - (spineLength - 1);
				std::cout << "UC: " << localId << std::endl;

				double parametricCoord = localId / (double)this->NumberOfRadialQuads;

				// Translation comes from the first point.
				input->GetPoint(spineIds[spineId][spineLength - 1], p0);

				// Axis of rotation comes from the centreline direction at the first point.
				input->GetPoint(spineIds[spineId][spineLength - 2], p1);
				vtkMath::Subtract(p1, p0, v0);

				// Get the last radius.
				radiiArray->GetTuple(spineIds[spineId][spineLength - 1], r);
				vtkMath::MultiplyScalar(r, -1.0);

				// Get the rotation axis.
				DoubleCross1(r, v0, r, c0);

				// Angle of rotation comes from the parametric coordinate along the arc.
				double angle = vtkMath::Pi() * parametricCoord;

				// Assemble local transform.
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->Translate(p0);
				transform->RotateWXYZ(vtkMath::DegreesFromRadians(angle), c0);

				// Transform the radius.
				transform->TransformPoint(r, point);

				// If not at the patch corner point, we need to insert a derivative.
				if(localId != 0)
				{
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
				// Inserting along the "right" patch edge.

				vtkIdType localId = ptId - this->NumberOfRadialQuads * 2  - (spineLength - 1);
				localId = std::fabs(localId - (spineLength - 1));
				std::cout << "RS: " << localId << std::endl;

				// Translation comes from the first point.
				input->GetPoint(spineIds[spineId][localId], p0);

				// Get the current radius.
				radiiArray->GetTuple(spineIds[spineId][localId], r);

				// Translate the radius.
				vtkMath::Add(p0, r, point);

				// If not at the patch corner point, we need to insert a derivative.
				if(localId != spineLength - 1)
				{
					// Get the centreline direction at the current point.
					input->GetPoint(spineIds[spineId][localId], p0);
					input->GetPoint(spineIds[spineId][localId + 1], p1);
					vtkMath::Subtract(p1, p0, v0);

					// The derivative is perpendicular to the radius vector and the local direction.
					vtkMath::Cross(v0, r, deriv);
					vtkMath::MultiplyScalar(deriv, -1.0);
					vtkMath::Normalize(deriv);

					double scaling = vtkMath::Norm(r) * yEdgeScaling;
					// Scale derivative vector magnitude.
					vtkMath::MultiplyScalar(deriv, scaling);
				}
			}

			vtkIdType id = patchPoints->InsertNextPoint(point);
			// Sanity check.
			assert(id == ptId);
			patchBoundary->GetPointIds()->InsertNextId(ptId);
			derivatives->InsertNextTuple(deriv);
		}
		patchBoundary->GetPointIds()->InsertNextId(0);

#if 1
		// Adjust the angles of derivatives for bifurcation segments.
		if(bifurcation)
		{
			radiiArray->GetTuple(spineIds[spineId][bifurcationPos], r);
			double scaling = vtkMath::Norm(r) * yEdgeScaling;

			double derivNm1[3];
			double derivNp1[3];
			double derivN[3];

			// For the "left" bifurcation point get the adjacent vectors.
			derivatives->GetTuple(rightBifurcationDerivId - 1, derivNm1);
			derivatives->GetTuple(rightBifurcationDerivId + 1, derivNp1);
			vtkMath::Add(derivNm1, derivNp1, derivN);
			vtkMath::Normalize(derivN);
			vtkMath::MultiplyScalar(derivN, scaling);
			derivatives->SetTuple(rightBifurcationDerivId, derivN);

			// Figure out the angles between the vectors at bifurcations and the adjacent vectors.
			double rightAngleNm1 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(derivNm1, derivN));
			double rightAngleNp1 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(derivNp1, derivN));

			// For the "right" bifurcation point get the adjacent vectors.
			derivatives->GetTuple(leftBifurcationDerivId - 1, derivNp1);
			derivatives->GetTuple(leftBifurcationDerivId + 1, derivNm1);
			vtkMath::Add(derivNm1, derivNp1, derivN);
			vtkMath::Normalize(derivN);
			vtkMath::MultiplyScalar(derivN, scaling);
			derivatives->SetTuple(leftBifurcationDerivId, derivN);

			// Figure out the angles between the vectors at bifurcations and the adjacent vectors.
			double leftAngleNm1 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(derivNm1, derivN));
			double leftAngleNp1 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(derivNp1, derivN));

			// TODO: This constant is to be examined closer.
			double rotationCoeff = 1.1;
#if 1
			// Traversing the "right" edge of the patch from bifurcation backwards.
			for(int ptId = bifurcationPos - 1, derivId = rightBifurcationDerivId - 1; ptId > 0; ptId--, derivId--)
			{
				// Fraction of the rotation angle.
				double currentFraction = ptId / (double)(bifurcationPos - 1);
				std::cout << "** LB; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

				// Rotate derivative around radius by fraction of the angle.
				radiiArray->GetTuple(spineIds[spineId][ptId], r);
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->RotateWXYZ(-rightAngleNm1 * rotationCoeff * currentFraction, r);

				std::cout << rightAngleNm1 * rotationCoeff * currentFraction << std::endl;

				double deriv0[3];
				derivatives->GetTuple(derivId, deriv0);
				double deriv1[3];
				transform->TransformPoint(deriv0, deriv1);
				derivatives->SetTuple(derivId, deriv1);
			}
			// Traversing the "left" edge of the patch from bifurcation backwards.
			for(int ptId = bifurcationPos - 1, derivId = leftBifurcationDerivId + 1; ptId > 0; ptId--, derivId++)
			{
				// Fraction of the rotation angle.
				double currentFraction = ptId / (double)(bifurcationPos - 1);
				std::cout << "++ LB; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

				// Rotate derivative around radius by fraction of the angle.
				radiiArray->GetTuple(spineIds[spineId][ptId], r);
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->RotateWXYZ(-leftAngleNm1 * rotationCoeff * currentFraction, r);

				std::cout << leftAngleNm1 * rotationCoeff * currentFraction << std::endl;

				double deriv0[3];
				derivatives->GetTuple(derivId, deriv0);
				double deriv1[3];
				transform->TransformPoint(deriv0, deriv1);
				derivatives->SetTuple(derivId, deriv1);
			}
			// Traversing the "right" edge of the patch from bifurcation forward.
			for(int ptId = bifurcationPos + 1, derivId = rightBifurcationDerivId + 1; ptId < spineLength - 1; ptId++, derivId++)
			{
				// Fraction of the rotation angle.
				double currentFraction = (spineLength - 1 - ptId) / (double)(spineLength - 1 - bifurcationPos);
				std::cout << "== LA; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

				// Rotate derivative around radius by fraction of the angle.
				radiiArray->GetTuple(spineIds[spineId][ptId], r);
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->RotateWXYZ(rightAngleNm1 * rotationCoeff * currentFraction, r);

				std::cout << rightAngleNm1 * rotationCoeff * currentFraction << std::endl;

				double deriv0[3];
				derivatives->GetTuple(derivId, deriv0);
				double deriv1[3];
				transform->TransformPoint(deriv0, deriv1);
				derivatives->SetTuple(derivId, deriv1);
			}
			// Traversing the "left" edge of the patch from bifurcation forward.
			for(int ptId = bifurcationPos + 1, derivId = leftBifurcationDerivId - 1; ptId < spineLength - 1; ptId++, derivId--)
			{
				// Fraction of the rotation angle.
				double currentFraction = (spineLength - 1 - ptId) / (double)(spineLength - 1 - bifurcationPos);
				std::cout << "-- RA; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

				// Rotate derivative around radius by fraction of the angle.
				radiiArray->GetTuple(spineIds[spineId][ptId], r);
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->RotateWXYZ(leftAngleNm1 * rotationCoeff * currentFraction, r);

				std::cout << leftAngleNm1 * rotationCoeff * currentFraction << std::endl;

				double deriv0[3];
				derivatives->GetTuple(derivId, deriv0);
				double deriv1[3];
				transform->TransformPoint(deriv0, deriv1);
				derivatives->SetTuple(derivId, deriv1);
			}
#endif
		}
#endif

		vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
		boundaries->InsertNextCell(patchBoundary);

		vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
		inputPatch->SetPoints(patchPoints);
		inputPatch->SetLines(boundaries);
		inputPatch->GetPointData()->SetVectors(derivatives);

		std::cout << inputPatch->GetNumberOfPoints() << " =?= " << numPtIds << std::endl;

		showPolyData(inputPatch, NULL, 0.1);

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

		vtkSmartPointer<vtkPolyData> outputPatch = vtkSmartPointer<vtkPolyData>::New();
		outputPatch->DeepCopy(patchFilter->GetOutput());
		vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();
		outputPoints->DeepCopy(outputPatch->GetPoints());

		vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
		structuredGrid->SetDimensions(this->NumberOfRadialQuads + 1, spineLength, 1);
		structuredGrid->SetPoints(outputPoints);

		//showPolyData(inputPatch, structuredGrid);
		outputGrids.push_back(structuredGrid);

		appendPolyDataFilter->AddInputData(inputPatch);

		//inputPatches.push_back(inputPatch);
	}

	appendPolyDataFilter->Update();
	showPolyData(appendPolyDataFilter->GetOutput(), NULL, 0.1);

	showGrids(outputGrids, input);

	// TODO: Merge output patches into one vtkPolyData object. vtkStructuredGridAppend filter from VTK nightly is not working correctly.

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
