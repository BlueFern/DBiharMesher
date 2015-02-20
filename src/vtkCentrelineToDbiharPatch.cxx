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

vtkStandardNewMacro(vtkCentrelineToDbiharPatch);


#if 0
#define SSTR( x ) dynamic_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()
#endif

vtkCentrelineToDbiharPatch::vtkCentrelineToDbiharPatch()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->cEdgeScaling = 4.0;
	this->yEdgeScaling = 6.0;
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

	int spineLength = pointIdList->GetNumberOfIds();

	// Find bifurcation point (if any).
	vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
	for (int ptPos = 0; ptPos < spineLength; ptPos++)
	{
		// TODO: Perhaps the loop doesn't need to be broken out of, just to make sure there is
		// a maximum of one bifurcation (if any).
		vtkIdType ptId = pointIdList->GetId(ptPos);
		input->GetPointCells(ptId, cells);

		std::cout << "### " << cells->GetNumberOfIds() << std::endl;

		if (cells->GetNumberOfIds() > 2) // Bifurcation point will involve (at least) 3 spines.
		{
			bifurcation = true;
			bifurcationPos = ptPos;
			//break;
		}
	}

	// Number of points in the patch boundary.
	int numPtIds = spineLength * 2 + this->NumberOfRadialQuads * 2 - 2;

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
	// The Diriclet boundary values are stored as the ordered set of points.
	// The Neumann boundary values are stored as vectors associated with each corresponding point.
	// What is called 'spineIds' in this case are just lists of point ids corresponding to segments of centrelines which are to be used for producing surface meshes.
	// For a non-bifurcating surface mesh there are two list of point ids, one per input patch for the Dbihar filter.
	// For a bifurcating surface mesh there are three lists of point ids, one per input patch for the Dbihar filter.
	// The code walks the given point id list starting at the first point. For this centreline point id an arch consisting of points obtained through radius rotation around this centreline point is inserted into the Dirichlet boundary condition point set.
	// Then the centreline point id list is traversed forward, where for each point id and the corresponding radius vector a point is added into the Dirichlet boundary condition point set.
	// When the last centreline point id is reached, an arch consisting of points obtained through radius rotation around this centreline point is inserted into the Dirichlet boundary condition point set.
	// Then the centreline point id list is traversed backcward, where for each point id and the corresponding radius vector a point is added into the Dirichlet boundary condition point set.
	// The Neumann boundary conditions are inserted for every point during the centreline point id list traversal. For the arches the direction of the derivatives is parallel to the direction of the centreline at the start and end points of the centreline segment.
	// For the other two edges of the patch the derivatives are oriented perpendicular to both the vector radius and the corresponding centreline edge. The magnitudes of the derivative vectors are proportional to the radius of the centreline at the given point.
	// While the code is in the debugg stage, the patches with the initialised boundary conditions are shown for visual verification.

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
				vtkMath::MultiplyScalar(deriv, cEdgeScaling);
			}
		}
		else if(ptId < this->NumberOfRadialQuads + spineLength - 1)
		{
			// Inserting along the "left" patch edge.

			vtkIdType localId = ptId - this->NumberOfRadialQuads;
			std::cout << "LS: " << localId << std::endl;

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
			input->GetPoint(pointIdList->GetId(spineLength - 1), p0);

			// Axis of rotation comes from the centreline direction at the first point.
			input->GetPoint(pointIdList->GetId(spineLength - 2), p1);
			vtkMath::Subtract(p1, p0, v0);

			// Get the last radius.
			radiiArray->GetTuple(pointIdList->GetId(spineLength - 1), r);
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
			input->GetPoint(pointIdList->GetId(localId), p0);

			// Get the current radius.
			radiiArray->GetTuple(pointIdList->GetId(localId), r);

			// Translate the radius.
			vtkMath::Add(p0, r, point);

			// If not at the patch corner point, we need to insert a derivative.
			if(localId != spineLength - 1)
			{
				// Get the centreline direction at the current point.
				input->GetPoint(pointIdList->GetId(localId), p0);
				input->GetPoint(pointIdList->GetId(localId + 1), p1);
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

	// Adjust the angles of derivatives for bifurcation segments.
	if (bifurcation)
	{
		radiiArray->GetTuple(pointIdList->GetId(bifurcationPos), r);
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
		double rightAngleNm1 = vtkMath::DegreesFromRadians(vtkDbiharStatic::AngleBetweenVectors(derivNm1, derivN));
		double rightAngleNp1 = vtkMath::DegreesFromRadians(vtkDbiharStatic::AngleBetweenVectors(derivNp1, derivN));

		// For the "right" bifurcation point get the adjacent vectors.
		derivatives->GetTuple(leftBifurcationDerivId - 1, derivNp1);
		derivatives->GetTuple(leftBifurcationDerivId + 1, derivNm1);
		vtkMath::Add(derivNm1, derivNp1, derivN);
		vtkMath::Normalize(derivN);
		vtkMath::MultiplyScalar(derivN, scaling);
		derivatives->SetTuple(leftBifurcationDerivId, derivN);

		// Figure out the angles between the vectors at bifurcations and the adjacent vectors.
		double leftAngleNm1 = vtkMath::DegreesFromRadians(vtkDbiharStatic::AngleBetweenVectors(derivNm1, derivN));
		double leftAngleNp1 = vtkMath::DegreesFromRadians(vtkDbiharStatic::AngleBetweenVectors(derivNp1, derivN));

		// TODO: This constant is to be examined closer.
		double rotationCoeff = 0.1;

		// Traversing the "right" edge of the patch from bifurcation backwards.
		for (int ptId = bifurcationPos - 1, derivId = rightBifurcationDerivId - 1; ptId > 0; ptId--, derivId--)
		{
			// Fraction of the rotation angle.
			double currentFraction = ptId / (double)(bifurcationPos - 1);
			std::cout << "** LB; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(pointIdList->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			double angle = -rightAngleNm1 * rotationCoeff * currentFraction;
			transform->RotateWXYZ(angle, r);
			std::cout << angle << std::endl;

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			double norm = vtkMath::Norm(deriv1);
			double norm1 = norm / cos(vtkMath::RadiansFromDegrees(std::fabs(angle)));
			double ratio = norm1/norm;
			vtkMath::MultiplyScalar(deriv1, ratio);
			std::cout << norm << ", " << norm1 << ", " << ratio << std::endl;

			derivatives->SetTuple(derivId, deriv1);
		}

		// Traversing the "left" edge of the patch from bifurcation backwards.
		for(int ptId = bifurcationPos - 1, derivId = leftBifurcationDerivId + 1; ptId > 0; ptId--, derivId++)
		{
			// Fraction of the rotation angle.
			double currentFraction = ptId / (double)(bifurcationPos - 1);
			std::cout << "++ LB; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(pointIdList->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			double angle = -leftAngleNm1 * rotationCoeff * currentFraction;
			transform->RotateWXYZ(angle, r);
			std::cout << angle << std::endl;

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			double norm = vtkMath::Norm(deriv1);
			double norm1 = norm / cos(vtkMath::RadiansFromDegrees(std::fabs(angle)));
			double ratio = norm1/norm;
			vtkMath::MultiplyScalar(deriv1, ratio);
			std::cout << norm << ", " << norm1 << ", " << ratio << std::endl;

			derivatives->SetTuple(derivId, deriv1);
		}

		// Traversing the "right" edge of the patch from bifurcation forward.
		for(int ptId = bifurcationPos + 1, derivId = rightBifurcationDerivId + 1; ptId < spineLength - 1; ptId++, derivId++)
		{
			// Fraction of the rotation angle.
			double currentFraction = (spineLength - 1 - ptId) / (double)(spineLength - 1 - bifurcationPos);
			std::cout << "== LA; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(pointIdList->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			double angle = rightAngleNm1 * rotationCoeff * currentFraction;
			transform->RotateWXYZ(angle, r);
			std::cout << angle << std::endl;

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			double norm = vtkMath::Norm(deriv1);
			double norm1 = norm / cos(vtkMath::RadiansFromDegrees(std::fabs(angle)));
			double ratio = norm1/norm;
			vtkMath::MultiplyScalar(deriv1, ratio);
			std::cout << norm << ", " << norm1 << ", " << ratio << std::endl;

			derivatives->SetTuple(derivId, deriv1);
		}

		// Traversing the "left" edge of the patch from bifurcation forward.
		for (int ptId = bifurcationPos + 1, derivId = leftBifurcationDerivId - 1; ptId < spineLength - 1; ptId++, derivId--)
		{
			// Fraction of the rotation angle.
			double currentFraction = (spineLength - 1 - ptId) / (double)(spineLength - 1 - bifurcationPos);
			std::cout << "-- RA; ptId: " << ptId << ", currentFraction: " << currentFraction << std::endl;

			// Rotate derivative around radius by fraction of the angle.
			radiiArray->GetTuple(pointIdList->GetId(ptId), r);
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			double angle = leftAngleNm1 * rotationCoeff * currentFraction;
			transform->RotateWXYZ(angle, r);
			std::cout << angle << std::endl;

			double deriv0[3];
			derivatives->GetTuple(derivId, deriv0);
			double deriv1[3];
			transform->TransformPoint(deriv0, deriv1);

			double norm = vtkMath::Norm(deriv1);
			double norm1 = norm / cos(vtkMath::RadiansFromDegrees(std::fabs(angle)));
			double ratio = norm1/norm;
			vtkMath::MultiplyScalar(deriv1, ratio);
			std::cout << norm << ", " << norm1 << ", " << ratio << std::endl;

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

	showPolyData(inputPatch, NULL, 0.1);

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	// Set the bounds of the UV space.
	patchFilter->SetA(0.0);
	patchFilter->SetB(1.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(1.0);
	// Set the boundary conditions.
	patchFilter->SetMQuads(this->NumberOfRadialQuads);
	patchFilter->SetNQuads(spineLength - 1);
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
