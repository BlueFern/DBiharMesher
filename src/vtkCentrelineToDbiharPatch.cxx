#include <vector>
#include <string>
#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkCallbackCommand.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkPointData.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkTransform.h>
#include <vtkIdList.h>
#include <vtkMath.h>

#include "vtkDbiharStatic.h"
#include "vtkCentrelineToDbiharPatch.h"
#include "vtkDbiharPatchFilter.h"

double expBase = 1.6;
double inputMap[] = {60, 7};
double fractionScale = 40;

double scaleMap(double input, double *inputMap)
{
	double input1 = inputMap[1] / inputMap[0] * std::abs(input);
	double expOut1 = std::pow(expBase, input1) - 1;
	double exp0ut2 = expOut1 / std::pow(expBase, inputMap[1]);
	return exp0ut2;
}

vtkStandardNewMacro(vtkCentrelineToDbiharPatch);

// TODO: Move this macro to a more appropriate place.
#if 1
#define SSTR( x ) dynamic_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()
#endif

vtkCentrelineToDbiharPatch::vtkCentrelineToDbiharPatch()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->NumberOfRadialQuads = 16;
	this->SpineId = 0;
	this->ArchDerivScale = 3;
	this->EdgeDerivScale = 4;

	this->ShowProgress = false;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkCentrelineToDbiharPatch::RequestData(vtkInformation *vtkNotUsed(request),
										vtkInformationVector **inputVector,
										vtkInformationVector *outputVector)
{
	vtkPolyData *input = vtkPolyData::GetData(inputVector[0],0);
	input->BuildLinks();

	vtkPolyData *output = vtkPolyData::GetData(outputVector,0);

	// TODO: Verify segment id and number of radial quads is set.

	// Find the spine indicated by the SpineId.
	vtkSmartPointer<vtkIdList> spineIds = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkCellArray> allSpines = input->GetLines();
	allSpines->InitTraversal();
	vtkIdType searchId = 0;
	while(true)
	{
		allSpines->GetNextCell(spineIds);
		if(searchId == this->SpineId)
		{
			break;
		}
		else
		{
			searchId++;
		}
	}
	const int spineSize = spineIds->GetNumberOfIds();

	// Get end points and bifurcations.
	vtkSmartPointer<vtkCellArray> verts = input->GetVerts();
	verts->InitTraversal();

	// Get end points for later.
	// WARNING: the assumption here is that end points always come first.
	vtkSmartPointer<vtkIdList> endPointIds = vtkSmartPointer<vtkIdList>::New();
	verts->GetNextCell(endPointIds);

	// Get bifurcation points.
	vtkSmartPointer<vtkIdList> bifurcationIds = vtkSmartPointer<vtkIdList>::New();
	verts->GetNextCell(bifurcationIds);

	// Determine if this is a bifurcation.
	// Remember the position if this is the case.
	// TODO: Check the spine contains only one bifurcation.
	bool bifurcation = false;
	vtkIdType bifurcationPos = -1;
	for(vtkIdType bifPos = 0; bifPos < bifurcationIds->GetNumberOfIds(); bifPos++)
	{
		bifurcationPos = spineIds->IsId(bifurcationIds->GetId(bifPos));
		if(bifurcationPos != -1)
		{
			bifurcation = true;
			break;
		}
	}

	// Points and lines.
	vtkSmartPointer<vtkPoints> patchPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> patchBoundary = vtkSmartPointer<vtkPolyLine>::New();

	// Derivatives.
	vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
	derivatives->SetName(vtkDbiharStatic::DERIV_ARR_NAME);
	derivatives->SetNumberOfComponents(3);

	// Number of points in the patch boundary.
	int numPtIds = spineSize * 2 + this->NumberOfRadialQuads * 2 - 2;

	vtkIdType rightBifurcationDerivId = -1;
	vtkIdType leftBifurcationDerivId = -1;

	// Every 10%, for progress function.
	int total = numPtIds / 10 ;
	int stage = 1;
	if (bifurcation)
	{
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

	// WARNING: There is a lot of code duplication here. This needs to be refactored as follows:
	// A single derivative vector is calculated for the lower arc and inserted to all arc points. Normalised.
	// A single derivative vector is calculated for the upper arc and inserted to all arc points. Normalised.
	// Traverse the spine and insert derivatives at each point in to a temporary container. Normalised.
	// Interpolate the angles in such a way that the angles between the derivative and the adjacent edges are equal.
	// Scale the derivatives in the temporary container.
	// Insert the derivatives from the temporary container into the right and left edges.
	// TODO: The above.

	vtkSmartPointer<vtkDoubleArray> radiiArray = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetVectors(vtkDbiharStatic::RADII_VECTORS_ARR_NAME));

	for(vtkIdType ptId = 0, spinePtId = 0; ptId < numPtIds; ptId++)
	{
		// Derivative is declared inside the loop to make sure it is (0,0,0) at the start of every iteration.
		double deriv[3] = {0.0};

		if(ptId < this->NumberOfRadialQuads)
		{
			// Inserting the first (lower edge) arc.

			vtkIdType localId = ptId;

			double parametricCoord = localId / (double)this->NumberOfRadialQuads;

			// Translation comes from the first point.
			input->GetPoint(spineIds->GetId(0), p0);

			// Axis of rotation comes from the centreline direction at the first point.
			input->GetPoint(spineIds->GetId(1), p1);
			vtkMath::Subtract(p1, p0, v0);

			// Get the first radius.
			radiiArray->GetTuple(spineIds->GetId(0), r);

			// Get the rotation axis.
			vtkDbiharStatic::DoubleCross(r, v0, r, c0);

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
				double scaling = vtkMath::Norm(r) * ArchDerivScale;
				vtkMath::MultiplyScalar(deriv, scaling);
			}
		}
		else if(ptId < this->NumberOfRadialQuads + spineSize - 1)
		{
			// Inserting along the "left" patch edge.

			vtkIdType localId = ptId - this->NumberOfRadialQuads;

			// Translation comes from the first point.
			input->GetPoint(spineIds->GetId(localId), p0);

			// Get the current radius.
			radiiArray->GetTuple(spineIds->GetId(localId), r);

			// Flip the radius.
			vtkMath::MultiplyScalar(r, -1.0);

			// Translate the radius.
			vtkMath::Add(p0, r, point);

			// If not at the patch corner point, we need to insert a derivative.
			if(localId != 0)
			{
				// Get the centreline direction at the current point.
				input->GetPoint(spineIds->GetId(localId), p0);
				input->GetPoint(spineIds->GetId(localId + 1), p1);
				vtkMath::Subtract(p1, p0, v0);

				// The derivative is perpendicular to the radius vector and the local direction.
				vtkMath::Cross(v0, r, deriv);
				vtkMath::Normalize(deriv);

				double scaling = vtkMath::Norm(r) * EdgeDerivScale;
				// Scale derivative vector magnitude.
				vtkMath::MultiplyScalar(deriv, scaling);
			}
		}
		else if(ptId < this->NumberOfRadialQuads * 2 + spineSize - 1)
		{
			// Inserting the second (upper edge) arc.

			vtkIdType localId = ptId - this->NumberOfRadialQuads - (spineSize - 1);

			double parametricCoord = localId / (double)this->NumberOfRadialQuads;

			// Translation comes from the first point.
			input->GetPoint(spineIds->GetId(spineSize - 1), p0);

			// Axis of rotation comes from the centreline direction at the first point.
			input->GetPoint(spineIds->GetId(spineSize - 2), p1);
			vtkMath::Subtract(p1, p0, v0);

			// Get the last radius.
			radiiArray->GetTuple(spineIds->GetId(spineSize - 1), r);
			vtkMath::MultiplyScalar(r, -1.0);

			// Get the rotation axis.
			vtkDbiharStatic::DoubleCross(r, v0, r, c0);

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
				double scaling = vtkMath::Norm(r) * ArchDerivScale;
				vtkMath::MultiplyScalar(deriv, scaling);
			}
		}
		else
		{
			// Inserting along the "right" patch edge.

			vtkIdType localId = ptId - this->NumberOfRadialQuads * 2  - (spineSize - 1);
			localId = std::fabs(localId - (spineSize - 1));

			// Translation comes from the first point.
			input->GetPoint(spineIds->GetId(localId), p0);

			// Get the current radius.
			radiiArray->GetTuple(spineIds->GetId(localId), r);

			// Translate the radius.
			vtkMath::Add(p0, r, point);

			// If not at the patch corner point, we need to insert a derivative.
			if(localId != spineSize - 1)
			{
				// Get the centreline direction at the current point.
				input->GetPoint(spineIds->GetId(localId), p0);
				input->GetPoint(spineIds->GetId(localId + 1), p1);
				vtkMath::Subtract(p1, p0, v0);

				// The derivative is perpendicular to the radius vector and the local direction.
				vtkMath::Cross(v0, r, deriv);
				vtkMath::MultiplyScalar(deriv, -1.0);
				vtkMath::Normalize(deriv);

				double scaling = vtkMath::Norm(r) * EdgeDerivScale;
				// Scale derivative vector magnitude.
				vtkMath::MultiplyScalar(deriv, scaling);
			}
		}

		vtkIdType id = patchPoints->InsertNextPoint(point);

		// Sanity check.
		assert(id == ptId);

		patchBoundary->GetPointIds()->InsertNextId(ptId);
		derivatives->InsertNextTuple(deriv);

		if ((ptId + 1) % total == 0)
		{
			// Over 12 calculated as 10% + 2 as we start with stage at 1 and have a later progress function call.
			this->UpdateProgress(static_cast<double>(stage++) / static_cast<double>(12));
		}
	}
	patchBoundary->GetPointIds()->InsertNextId(0);

	vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
	boundaries->InsertNextCell(patchBoundary);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(patchPoints);
	inputPatch->SetLines(boundaries);
	inputPatch->GetPointData()->SetVectors(derivatives);

	// vtkDbiharStatic::ShowPolyData(inputPatch);

	// Adjust the angles of derivatives for bifurcation segments.
	if(bifurcation)
	{

		// TODO: Perhaps this constant needs to be a parameter?
		const double rotationCoeff = 0.1;

		double derivNm1[3];
		double derivNp1[3];
		double derivN[3];

		// For the "left" bifurcation point get the adjacent vectors.
		derivatives->GetTuple(rightBifurcationDerivId - 1, derivNm1);
		derivatives->GetTuple(rightBifurcationDerivId + 1, derivNp1);
		// Take the mean of those vectors.
		vtkMath::Add(derivNm1, derivNp1, derivN);
		vtkMath::Normalize(derivN);

		// Figure out the angles between the vectors at bifurcations and the adjacent vectors.
		double rightAngleNm1 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(derivNm1, derivN));
		double rightAngleNp1 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(derivNp1, derivN));

		vtkMath::MultiplyScalar(derivN, vtkMath::Norm(derivatives->GetTuple(rightBifurcationDerivId)) * 1.6); // * scalingR); // * this->EdgeDerivScale);
		derivatives->SetTuple(rightBifurcationDerivId, derivN);

		// For the "right" bifurcation point get the adjacent vectors.
		derivatives->GetTuple(leftBifurcationDerivId - 1, derivNp1);
		derivatives->GetTuple(leftBifurcationDerivId + 1, derivNm1);
		// Take the mean of those vectors.
		vtkMath::Add(derivNm1, derivNp1, derivN);
		vtkMath::Normalize(derivN);

		// Figure out the angles between the vectors at bifurcations and the adjacent vectors.
		double leftAngleNm1 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(derivNm1, derivN));
		double leftAngleNp1 = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(derivNp1, derivN));

		vtkMath::MultiplyScalar(derivN, vtkMath::Norm(derivatives->GetTuple(leftBifurcationDerivId)) * 1.6); // * scalingL); // * this->EdgeDerivScale);
		derivatives->SetTuple(leftBifurcationDerivId, derivN);

		// vtkDbiharStatic::ShowPolyData(inputPatch);

		// Derivatives adjustment rotation axis.
		radiiArray->GetTuple(spineIds->GetId(bifurcationPos), r);

		double m = 4.0 / 60.0;

		// Traversing the "right" edge of the patch from bifurcation backwards.
		for (int ptId = bifurcationPos - 1, derivId = rightBifurcationDerivId - 1; ptId > 0; ptId--, derivId--)
		{
			// Fraction of the rotation angle.
			double currentFraction = ptId / (double)(bifurcationPos - 1);
			std::cout << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(spineIds->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			double angle = -rightAngleNm1 * currentFraction;
			angle += -rightAngleNm1 * m * currentFraction;
			transform->RotateWXYZ(angle, r);

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			vtkMath::Normalize(deriv1);
			vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(angle, inputMap));
			//std::cout << "\t" << vtkMath::Norm(deriv0) << ", " << vtkMath::Norm(deriv0) * scaleMap(angle, inputMap) << ", " << angle << std::endl;

			//vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(currentFraction * fractionScale, inputMap));
			//std::cout << "\t" << vtkMath::Norm(deriv0) << ", " << vtkMath::Norm(deriv0) * scaleMap(currentFraction * fractionScale, inputMap) << ", " << currentFraction << std::endl;

			derivatives->SetTuple(derivId, deriv1);
		}

		std::cout << std::endl;


		// Traversing the "left" edge of the patch from bifurcation backwards.
		for(int ptId = bifurcationPos - 1, derivId = leftBifurcationDerivId + 1; ptId > 0; ptId--, derivId++)
		{
			// Fraction of the rotation angle.
			double currentFraction = ptId / (double)(bifurcationPos - 1);
			std::cout << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(spineIds->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			double angle = -leftAngleNm1 * currentFraction;
			angle += -leftAngleNm1 * m * currentFraction;
			transform->RotateWXYZ(angle, r);

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			vtkMath::Normalize(deriv1);
			vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(angle, inputMap));
			//std::cout << "\t" << vtkMath::Norm(deriv0) << ", " << vtkMath::Norm(deriv0) * scaleMap(angle, inputMap) << ", " << angle << std::endl;

			//vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(currentFraction * fractionScale, inputMap));
			//std::cout << "\t" << vtkMath::Norm(deriv0) << ", " << vtkMath::Norm(deriv0) * scaleMap(currentFraction * fractionScale, inputMap) << ", " << currentFraction << std::endl;

			derivatives->SetTuple(derivId, deriv1);
		}

		std::cout << std::endl;

		// Traversing the "right" edge of the patch from bifurcation forward.
		for(int ptId = bifurcationPos + 1, derivId = rightBifurcationDerivId + 1; ptId < spineSize - 1; ptId++, derivId++)
		{
			// Fraction of the rotation angle.
			double currentFraction = (spineSize - ptId) / (double)(spineSize - 1 - bifurcationPos);
			std::cout << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(spineIds->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			double angle = rightAngleNm1 * currentFraction;
			angle += rightAngleNm1 * m * currentFraction;
			transform->RotateWXYZ(angle, r);

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			vtkMath::Normalize(deriv1);
			vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(angle, inputMap));
			//std::cout << "\t" << vtkMath::Norm(deriv0) << ", " << vtkMath::Norm(deriv0) * scaleMap(angle, inputMap) << ", " << angle << std::endl;

			//vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(currentFraction * fractionScale, inputMap));
			//std::cout << "\t" << vtkMath::Norm(deriv0) << ", " << vtkMath::Norm(deriv0) * scaleMap(currentFraction * fractionScale, inputMap) << ", " << currentFraction << std::endl;

			derivatives->SetTuple(derivId, deriv1);
		}

		std::cout << std::endl;

		// Traversing the "left" edge of the patch from bifurcation forward.
		for (int ptId = bifurcationPos + 1, derivId = leftBifurcationDerivId - 1; ptId < spineSize - 1; ptId++, derivId--)
		{
			// Fraction of the rotation angle.
			double currentFraction = (spineSize - ptId) / (double)(spineSize - 1 - bifurcationPos);
			std::cout << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(spineIds->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();

			double angle = leftAngleNm1 * currentFraction;
			angle += leftAngleNm1 * m * currentFraction;
			transform->RotateWXYZ(angle, r);

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			vtkMath::Normalize(deriv1);
			vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(angle, inputMap));
			//std::cout << "\t" << vtkMath::Norm(deriv0) << ", " << vtkMath::Norm(deriv0) * scaleMap(angle, inputMap) << ", " << angle << std::endl;

			//vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(currentFraction * fractionScale, inputMap));
			//std::cout << "\t" << vtkMath::Norm(deriv0) << ", " << vtkMath::Norm(deriv0) * scaleMap(currentFraction * fractionScale, inputMap) << ", " << currentFraction << std::endl;

			derivatives->SetTuple(derivId, deriv1);
		}
	}

	this->UpdateProgress(static_cast<double>(stage++) / static_cast<double>(12));

	// Sanity check.
	assert(inputPatch->GetNumberOfPoints() == numPtIds);

	// vtkDbiharStatic::ShowPolyData(inputPatch, 1.0);

	// writePolyData(inputPatch, SSTR("tmpPatch" << this->SpineId << ".vtp"));

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	double lowerArcLength = vtkMath::Pi() * vtkMath::Norm(radiiArray->GetTuple(spineIds->GetId(0)));
	double spineLength = 0;
	for(int pos = 0; pos < spineSize - 2; pos++)
	{
		input->GetPoint(spineIds->GetId(pos), p0);
		input->GetPoint(spineIds->GetId(pos + 1), p1);
		spineLength += sqrt(vtkMath::Distance2BetweenPoints(p0, p1));
	}

	// Set the bounds of the UV space.
	patchFilter->SetA(0.0);
	patchFilter->SetB(1.0);
	patchFilter->SetC(0.0);
	// patchFilter->SetD(1.0);
	patchFilter->SetD(spineLength / lowerArcLength);
	// Set the boundary conditions.
	patchFilter->SetMQuads(this->NumberOfRadialQuads);
	patchFilter->SetNQuads(spineSize - 1);
	// Set solution method.
	patchFilter->SetIFlag(2);

	patchFilter->SetInputData(inputPatch);
	patchFilter->Update();

	output->DeepCopy(patchFilter->GetOutput());
	return 1;
}

void vtkCentrelineToDbiharPatch::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkCentrelineToDbiharPatch* filter = static_cast<vtkCentrelineToDbiharPatch *>(caller);
	if(filter->ShowProgress)
	{
		cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
	}
}

void vtkCentrelineToDbiharPatch::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);

	os << indent << "SpineId: " << this->SpineId << "\n";
	os << indent << "ArchDerivScale: " << this->ArchDerivScale << "\n";
	os << indent << "EdgeDerivScale: " << this->EdgeDerivScale << "\n";
	os << indent << "NumberOfRadialQuads: " << this->NumberOfRadialQuads << "\n";
}
