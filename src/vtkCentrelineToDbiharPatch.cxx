#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkCallbackCommand.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkDataObject.h>
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
#include "showPolyData.h"

double expBase = 1.6;

double inputMap[] = {60, 7.0};

double scaleMap(double input, double *inputMap)
{
	double input1 = inputMap[1] / inputMap[0] * std::abs(input);

	double expOut1 = std::pow(expBase, input1) - 1;

	double exp0ut2 = expOut1 / std::pow(expBase, inputMap[1]);

	return exp0ut2;
}

#if 0
double abc0[] = {0.000493, 0, 0};

double quadraticFunction(double x, double *abc)
{
	return abc[0] * x * x + abc[1] * x + abc[2];
}
#endif

vtkStandardNewMacro(vtkCentrelineToDbiharPatch);

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
}

int vtkCentrelineToDbiharPatch::RequestData(vtkInformation *vtkNotUsed(request),
										vtkInformationVector **inputVector,
										vtkInformationVector *outputVector)
{
	vtkPolyData *input = vtkPolyData::GetData(inputVector[0],0);
	input->BuildLinks();

	vtkPolyData *output = vtkPolyData::GetData(outputVector,0);

	// TODO: Verify segment id and number of radial quads is set.

	bool bifurcation = false;
	vtkIdType bifurcationPos = -1;
	std::vector<vtkSmartPointer<vtkStructuredGrid> > outputGrids;

	vtkSmartPointer<vtkAppendPolyData> appendPolyDataFilter = vtkSmartPointer<vtkAppendPolyData>::New();
	vtkSmartPointer<vtkAppendPolyData> appendGridFilter = vtkSmartPointer<vtkAppendPolyData>::New();

	vtkSmartPointer<vtkDoubleArray> radiiArray =
			vtkDoubleArray::SafeDownCast(input->GetPointData()->GetVectors(vtkDbiharStatic::RADII_VECTORS_ARR_NAME));

	// Points and lines.
	vtkSmartPointer<vtkPoints> patchPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> patchBoundary = vtkSmartPointer<vtkPolyLine>::New();

	// Derivatives.
	vtkSmartPointer<vtkDoubleArray> derivatives = vtkSmartPointer<vtkDoubleArray>::New();
	derivatives->SetName(vtkDbiharStatic::DERIV_ARR_NAME);
	derivatives->SetNumberOfComponents(3);

	vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
	pointIdList = input->GetCell(this->SpineId)->GetPointIds();

	const int spineSize = pointIdList->GetNumberOfIds();

	// Find bifurcation point (if any).
	vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
	for (int ptPos = 0; ptPos < spineSize; ptPos++)
	{
		// TODO: Perhaps the loop doesn't need to be broken out of, just to make sure there is
		// a maximum of one bifurcation (if any).
		vtkIdType ptId = pointIdList->GetId(ptPos);
		input->GetPointCells(ptId, cells);

		// TODO: This is a temporary work around here. Bifurcation points need to be stored as a collection of vertices.
		if (cells->GetNumberOfIds() == 3) // Bifurcation point will involve (at least) 3 spines.
		{
			bifurcation = true;
			bifurcationPos = ptPos;
			// TODO: Need to go through all points in this spine and make sure there is not more than one bifurcation.
			// Throw an error (vtkErrorMacro) if this is not the case.
			// break;
		}
	}

	// Number of points in the patch boundary.
	int numPtIds = spineSize * 2 + this->NumberOfRadialQuads * 2 - 2;

	vtkIdType rightBifurcationDerivId = -1;
	vtkIdType leftBifurcationDerivId = -1;

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

	// TODO: This content should appear in the filter documentation.

	// Set Dirichlet and Neumann boundary conditions. Both boundary conditions are stored in a vtkPolyData object.
	// The Dirichlet boundary values are stored as the ordered set of points.
	// The Neumann boundary values are stored as vectors associated with each corresponding point.
	// What is called 'spineIds' in this case are just lists of point ids corresponding to segments of centrelines which are to be used for producing surface meshes.
	// For a non-bifurcating surface mesh there are two list of point ids, one per input patch for the Dbihar filter.
	// For a bifurcating surface mesh there are three lists of point ids, one per input patch for the Dbihar filter.
	// The code walks the given point id list starting at the first point. For this centreline point id an arch consisting of points obtained through radius rotation around this centreline point is inserted into the Dirichlet boundary condition point set.
	// Then the centreline point id list is traversed forward, where for each point id and the corresponding radius vector a point is added into the Dirichlet boundary condition point set.
	// When the last centreline point id is reached, an arch consisting of points obtained through radius rotation around this centreline point is inserted into the Dirichlet boundary condition point set.
	// Then the centreline point id list is traversed backward, where for each point id and the corresponding radius vector a point is added into the Dirichlet boundary condition point set.
	// The Neumann boundary conditions are inserted for every point during the centreline point id list traversal. For the arches the direction of the derivatives is parallel to the direction of the centreline at the start and end points of the centreline segment.
	// For the other two edges of the patch the derivatives are oriented perpendicular to both the vector radius and the corresponding centreline edge. The magnitudes of the derivative vectors are proportional to the radius of the centreline at the given point.
	// While the code is in the debug stage, the patches with the initialised boundary conditions are shown for visual verification.

	// TODO: Interpolate the derivatives (similar to what we do with radii) between adjacent points. Alternatively, we should try
	// calculating the derivatives only for the centreline and translating them to both edges. They would still have to be interpolated.

	// WARNING: There is a lot of code duplication here. This needs to be refactored as follows:
	// A single derivative vector is calculated for the lower arc and inserted to all arc points. Normalised.
	// A single derivative vector is calculated for the upper arc and inserted to all arc points. Normalised.
	// Traverse the spine and insert derivatives at each point in to a temporary container. Normalised.
	// Interpolate the angles in such a way that the angles between the derivative and the adjacent edges are equal.
	// Scale the derivatives in the temporary container.
	// Insert the derivatives from the temporary container into the right and left edges.

	for(vtkIdType ptId = 0, spinePtId = 0; ptId < numPtIds; ptId++)
	{
		// Derivative is declared inside the loop to make sure it is (0,0,0) at the start of every iteration.
		double deriv[3] = {0.0};

		if(ptId < this->NumberOfRadialQuads)
		{
			// Inserting the first (lower edge) arc.

			vtkIdType localId = ptId;
			//// std::cout << "LC: " << localId << std::endl;

			double parametricCoord = localId / (double)this->NumberOfRadialQuads;

			// Translation comes from the first point.
			input->GetPoint(pointIdList->GetId(0), p0);

			// Axis of rotation comes from the centreline direction at the first point.
			input->GetPoint(pointIdList->GetId(1), p1);
			vtkMath::Subtract(p1, p0, v0);

			// Get the first radius.
			radiiArray->GetTuple(pointIdList->GetId(0), r);

			// Get the rotation axis.
			vtkDbiharStatic::DoubleCross1(r, v0, r, c0);

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
				double scaling = vtkMath::Norm(r) * ArchDerivScale;
				// double scaling = vtkMath::Norm(radiiArray->GetTuple(0)) * cEdgeScaling;
				vtkMath::MultiplyScalar(deriv, scaling);
			}
		}
		else if(ptId < this->NumberOfRadialQuads + spineSize - 1)
		{
			// Inserting along the "left" patch edge.

			vtkIdType localId = ptId - this->NumberOfRadialQuads;
			//// std::cout << "LS: " << localId << std::endl;

			// Translation comes from the first point.
			input->GetPoint(pointIdList->GetId(localId), p0);

			// Get the current radius.
			radiiArray->GetTuple(pointIdList->GetId(localId), r);

			// Flip the radius.
			vtkMath::MultiplyScalar(r, -1.0);

			// Translate the radius.
			vtkMath::Add(p0, r, point);

			// If not at the patch corner point, we need to insert a derivative.
			if(localId != 0)
			{
				// Get the centreline direction at the current point.
				input->GetPoint(pointIdList->GetId(localId), p0);
				input->GetPoint(pointIdList->GetId(localId + 1), p1);
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
			//// std::cout << "UC: " << localId << std::endl;

			double parametricCoord = localId / (double)this->NumberOfRadialQuads;

			// Translation comes from the first point.
			input->GetPoint(pointIdList->GetId(spineSize - 1), p0);

			// Axis of rotation comes from the centreline direction at the first point.
			input->GetPoint(pointIdList->GetId(spineSize - 2), p1);
			vtkMath::Subtract(p1, p0, v0);

			// Get the last radius.
			radiiArray->GetTuple(pointIdList->GetId(spineSize - 1), r);
			vtkMath::MultiplyScalar(r, -1.0);

			// Get the rotation axis.
			vtkDbiharStatic::DoubleCross1(r, v0, r, c0);

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
				double scaling = vtkMath::Norm(r) * ArchDerivScale;
				// double scaling = vtkMath::Norm(radiiArray->GetTuple(0)) * cEdgeScaling;
				vtkMath::MultiplyScalar(deriv, scaling);
			}
		}
		else
		{
			// Inserting along the "right" patch edge.

			vtkIdType localId = ptId - this->NumberOfRadialQuads * 2  - (spineSize - 1);
			localId = std::fabs(localId - (spineSize - 1));
			//// std::cout << "RS: " << localId << std::endl;

			// Translation comes from the first point.
			input->GetPoint(pointIdList->GetId(localId), p0);

			// Get the current radius.
			radiiArray->GetTuple(pointIdList->GetId(localId), r);

			// Translate the radius.
			vtkMath::Add(p0, r, point);

			// If not at the patch corner point, we need to insert a derivative.
			if(localId != spineSize - 1)
			{
				// Get the centreline direction at the current point.
				input->GetPoint(pointIdList->GetId(localId), p0);
				input->GetPoint(pointIdList->GetId(localId + 1), p1);
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
		//// std::cout << "*** Deriv norm: " << vtkMath::Norm(deriv) << std::endl;
	}
	patchBoundary->GetPointIds()->InsertNextId(0);

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
		// vtkMath::MultiplyScalar(derivN, vtkMath::Norm(derivatives->GetTuple(rightBifurcationDerivId)));
		// Store the scaled derivative vector.
		// derivatives->SetTuple(rightBifurcationDerivId, derivN);

		// Figure out the angles between the vectors at bifurcations and the adjacent vectors.
		double rightAngleNm1 = vtkMath::DegreesFromRadians(vtkDbiharStatic::AngleBetweenVectors(derivNm1, derivN));
		double rightAngleNp1 = vtkMath::DegreesFromRadians(vtkDbiharStatic::AngleBetweenVectors(derivNp1, derivN));

		//double tmpR = ((rightAngleNm1 + rightAngleNp1) * rotationCoeff) / 2.0;
		//double scalingR = quadraticFunction((rightAngleNm1 + rightAngleNp1) / 2.0, abc0);
		//double scalingR = scaleMap((rightAngleNm1 + rightAngleNp1) / 2.0, inputMap);
		vtkMath::MultiplyScalar(derivN, vtkMath::Norm(derivatives->GetTuple(rightBifurcationDerivId)) + vtkMath::Norm(derivatives->GetTuple(rightBifurcationDerivId))); // * scalingR); // * this->EdgeDerivScale);
		derivatives->SetTuple(rightBifurcationDerivId, derivN);

		// For the "right" bifurcation point get the adjacent vectors.
		derivatives->GetTuple(leftBifurcationDerivId - 1, derivNp1);
		derivatives->GetTuple(leftBifurcationDerivId + 1, derivNm1);
		// Take the mean of those vectors.
		vtkMath::Add(derivNm1, derivNp1, derivN);
		vtkMath::Normalize(derivN);
		// vtkMath::MultiplyScalar(derivN, vtkMath::Norm(derivatives->GetTuple(leftBifurcationDerivId)));
		// Store the scaled derivative vector.
		// derivatives->SetTuple(leftBifurcationDerivId, derivN);

		// Figure out the angles between the vectors at bifurcations and the adjacent vectors.
		double leftAngleNm1 = vtkMath::DegreesFromRadians(vtkDbiharStatic::AngleBetweenVectors(derivNm1, derivN));
		double leftAngleNp1 = vtkMath::DegreesFromRadians(vtkDbiharStatic::AngleBetweenVectors(derivNp1, derivN));

		//double tmpL = ((leftAngleNm1 + leftAngleNp1) * rotationCoeff) / 2.0;
		//double scalingL = quadraticFunction((leftAngleNm1 + leftAngleNp1) / 2.0, abc0);
		//double scalingL = scaleMap((leftAngleNm1 + leftAngleNp1) / 2.0, inputMap);
		vtkMath::MultiplyScalar(derivN, vtkMath::Norm(derivatives->GetTuple(leftBifurcationDerivId)) + vtkMath::Norm(derivatives->GetTuple(leftBifurcationDerivId))); // * scalingL); // * this->EdgeDerivScale);
		derivatives->SetTuple(leftBifurcationDerivId, derivN);

		//std::cout << "**************** tmpR: " << tmpR << ", tmpL: " << tmpL << std::endl;
		//std::cout << "******************** scalingR: " << scalingR << ", scalingL: " << scalingL << std::endl;

		std::cout << "*** " << this->SpineId << std::endl;
		std::cout << rightAngleNm1 << ", " << rightAngleNp1 << ", " << leftAngleNm1 << ", " << leftAngleNp1 << std::endl;

		radiiArray->GetTuple(pointIdList->GetId(bifurcationPos), r);

		double m = 4.0 / 60.0;

		// Traversing the "right" edge of the patch from bifurcation backwards.
		for (int ptId = bifurcationPos - 1, derivId = rightBifurcationDerivId - 1; ptId > 0; ptId--, derivId--)
		{
			// Fraction of the rotation angle.
			double currentFraction = ptId / (double)(bifurcationPos - 1);
			//// std::cout << "** LB; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(pointIdList->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			// double angle = -rightAngleNm1 * rotationCoeff * currentFraction;
			double angle = -rightAngleNm1 * currentFraction;
			angle += -rightAngleNm1 * m * currentFraction;
			transform->RotateWXYZ(angle, r);
			//// std::cout << angle << std::endl;

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			//double norm = vtkMath::Norm(deriv1);
			//double norm1 = norm / cos(vtkMath::RadiansFromDegrees(std::fabs(angle)));
			//double ratio = norm1/norm;
			//vtkMath::MultiplyScalar(deriv1, ratio);
			//// std::cout << norm << ", " << norm1 << ", " << ratio << std::endl;
			//std::cout << "ratio: " << ratio << ", angle: " << angle << std::endl;
			//std::cout << "angle: " << angle << std::endl;

			vtkMath::Normalize(deriv1);
			vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(angle, inputMap));

			derivatives->SetTuple(derivId, deriv1);
		}

		// Traversing the "left" edge of the patch from bifurcation backwards.
		for(int ptId = bifurcationPos - 1, derivId = leftBifurcationDerivId + 1; ptId > 0; ptId--, derivId++)
		{
			// Fraction of the rotation angle.
			double currentFraction = ptId / (double)(bifurcationPos - 1);
			//// std::cout << "++ LB; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(pointIdList->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			// double angle = -leftAngleNm1 * rotationCoeff * currentFraction;
			double angle = -leftAngleNm1 * currentFraction;
			angle += -leftAngleNm1 * m * currentFraction;
			transform->RotateWXYZ(angle, r);
			//// std::cout << angle << std::endl;

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			//double norm = vtkMath::Norm(deriv1);
			//double norm1 = norm / cos(vtkMath::RadiansFromDegrees(std::fabs(angle)));
			//double ratio = norm1/norm;
			//vtkMath::MultiplyScalar(deriv1, ratio);
			//// std::cout << norm << ", " << norm1 << ", " << ratio << std::endl;
			//std::cout << "ratio: " << ratio << ", angle: " << angle << std::endl;
			//std::cout << "angle: " << angle << std::endl;

			vtkMath::Normalize(deriv1);
			vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(angle, inputMap));

			derivatives->SetTuple(derivId, deriv1);
		}

		// Traversing the "right" edge of the patch from bifurcation forward.
		for(int ptId = bifurcationPos + 1, derivId = rightBifurcationDerivId + 1; ptId < spineSize - 1; ptId++, derivId++)
		{
			// Fraction of the rotation angle.
			double currentFraction = (spineSize - 1 - ptId) / (double)(spineSize - 1 - bifurcationPos);
			//// std::cout << "== LA; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(pointIdList->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			// double angle = rightAngleNm1 * rotationCoeff * currentFraction;
			double angle = rightAngleNm1 * currentFraction;
			angle += rightAngleNm1 * m * currentFraction;
			transform->RotateWXYZ(angle, r);
			//// std::cout << angle << std::endl;

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			//double norm = vtkMath::Norm(deriv1);
			//double norm1 = norm / cos(vtkMath::RadiansFromDegrees(std::fabs(angle)));
			//double ratio = norm1/norm;
			//vtkMath::MultiplyScalar(deriv1, ratio);
			//// std::cout << norm << ", " << norm1 << ", " << ratio << std::endl;
			//std::cout << "ratio: " << ratio << ", angle: " << angle << std::endl;
			//std::cout << "angle: " << angle << std::endl;

			vtkMath::Normalize(deriv1);
			vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(angle, inputMap));

			derivatives->SetTuple(derivId, deriv1);
		}

		// Traversing the "left" edge of the patch from bifurcation forward.
		for (int ptId = bifurcationPos + 1, derivId = leftBifurcationDerivId - 1; ptId < spineSize - 1; ptId++, derivId--)
		{
			// Fraction of the rotation angle.
			double currentFraction = (spineSize - 1 - ptId) / (double)(spineSize - 1 - bifurcationPos);
			//// std::cout << "-- RA; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(pointIdList->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			// double angle = leftAngleNm1 * rotationCoeff * currentFraction;
			double angle = leftAngleNm1 * currentFraction;
			angle += leftAngleNm1 * m * currentFraction;
			transform->RotateWXYZ(angle, r);
			//// std::cout << angle << std::endl;

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			//double norm = vtkMath::Norm(deriv1);
			//double norm1 = norm / cos(vtkMath::RadiansFromDegrees(std::fabs(angle)));
			//double ratio = norm1/norm;
			//vtkMath::MultiplyScalar(deriv1, ratio);
			//// std::cout << norm << ", " << norm1 << ", " << ratio << std::endl;
			//std::cout << "ratio: " << ratio << ", angle: " << angle << std::endl;
			//std::cout << "angle: " << angle << std::endl;

			vtkMath::Normalize(deriv1);
			vtkMath::MultiplyScalar(deriv1, vtkMath::Norm(deriv0) + vtkMath::Norm(deriv0) * scaleMap(angle, inputMap));

			derivatives->SetTuple(derivId, deriv1);
		}
	}

	vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
	boundaries->InsertNextCell(patchBoundary);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(patchPoints);
	inputPatch->SetLines(boundaries);
	inputPatch->GetPointData()->SetVectors(derivatives);

	// Sanity check.
	assert(inputPatch->GetNumberOfPoints() == numPtIds);

	// showPolyData(inputPatch, NULL, 1.0);

	// writePolyData(inputPatch, SSTR("tmpPatch" << this->SpineId << ".vtp"));

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	double lowerArcLength = vtkMath::Pi() * vtkMath::Norm(radiiArray->GetTuple(pointIdList->GetId(0)));
	double spineLength = 0;
	for(int pos = 0; pos < spineSize - 2; pos++)
	{
		input->GetPoint(pointIdList->GetId(pos), p0);
		input->GetPoint(pointIdList->GetId(pos + 1), p1);
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


/*


	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetDimensions(this->NumberOfRadialQuads + 1, spineLength, 1);
	structuredGrid->SetPoints(patchFilter->GetOutput()->GetPoints());

	vtkSmartPointer<vtkStructuredGridGeometryFilter> gridToPolyDataFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	gridToPolyDataFilter->SetInputData(structuredGrid);
	std::cout << this->NumberOfRadialQuads + 1 << ", " << spineLength << std::endl;
	std::cout << "points in grid: " << structuredGrid->GetNumberOfPoints() << std::endl;

	appendGridFilter->AddInputConnection(gridToPolyDataFilter->GetOutputPort());

	#if 0
	vtkSmartPointer<vtkXMLPolyDataWriter> pdw = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	pdw->SetInputData(patchFilter->GetOutput());
	pdw->SetFileName(SSTR("dbhPatch" << spineId << ".vtp").c_str());
	pdw->Write();

	vtkSmartPointer<vtkXMLStructuredGridWriter> writer = vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
	writer->SetFileName(SSTR("structuredGrid" << spineId << ".vts").c_str());
	writer->SetInputData(structuredGrid);
	writer->Update();
	#endif

	outputGrids.push_back(structuredGrid);

	appendPolyDataFilter->AddInputData(inputPatch);

	#if 0
	// TODO: Remove this.
	appendPolyDataFilter->Update();
	showPolyData(appendPolyDataFilter->GetOutput(), NULL, 0.1);

	appendGridFilter->Update();
	std::cout << "total points in grids: " << appendGridFilter->GetOutput()->GetNumberOfPoints() << std::endl;

	vtkSmartPointer<vtkCleanPolyData> polyDataCleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	polyDataCleaner->SetInputData(appendGridFilter->GetOutput());
	polyDataCleaner->Update();

	std::cout << "points after clean: " << polyDataCleaner->GetOutput()->GetNumberOfPoints() << std::endl;

	showGrids(outputGrids, input);
	#endif
	return 1;
*/
}

void vtkCentrelineToDbiharPatch::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	 // TODO: print class member values.
	os << indent << "\n";
}
