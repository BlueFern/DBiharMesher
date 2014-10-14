#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkIdList.h>
#include <vtkPolyLine.h>
#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkTransform.h>
#include <vtkStructuredGrid.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkMath.h>

#include "vtkDbiharPatchFilter.h"

#include "showPolyData.h"

double returnPoint[3] = {0.0};
void NegateVector(const double v[3], double nV[3])
{
	nV[0] = -v[0];
	nV[1] = -v[1];
	nV[2] = -v[2];
}
double* NegateVector(const double v[3])
{
	NegateVector(v, returnPoint);
	return returnPoint;
}

// This is taken from vtkMath class in VTK Nightly in October 2014.
double AngleBetweenVectors(const double v1[3], const double v2[3])
{
  double cross[3];
  vtkMath::Cross(v1, v2, cross);
  return atan2(vtkMath::Norm(cross), vtkMath::Dot(v1, v2));
}

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	int numPoints = 61;
	int numTrunkPoints = 31;
	int numBranchPoints = 30;

	double stepLen = 10.0;
	double radius = 40.0;

	vtkSmartPointer<vtkPoints> spinePoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> spine = vtkSmartPointer<vtkPolyLine>::New();
	vtkIdList *spineIds = spine->GetPointIds();

	double *radiiValues = new double[numPoints];
	std::fill_n(radiiValues, numPoints, radius);
	vtkSmartPointer<vtkDoubleArray> spineRadii = vtkSmartPointer<vtkDoubleArray>::New();
	spineRadii->SetArray(radiiValues, numPoints, 1);

	vtkSmartPointer<vtkCellArray> spineLines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTransform> spinePointTransform = vtkSmartPointer<vtkTransform>::New();

	double rootPoint[3] = {0.0, 0.0, 0.0};
	double inPoint[3]; memcpy(inPoint, rootPoint, sizeof(double) * 3); // Initialise to rootPoint.
	double stepVector[3] = {rootPoint[0], rootPoint[1] + stepLen, rootPoint[2]};

	for(int i = 0; i < numTrunkPoints; i++)
	{
		double outPoint1[3];
		spinePointTransform->TransformPoint(rootPoint, outPoint1);
		vtkIdType pId0 = spinePoints->InsertNextPoint(outPoint1);
		std::cout << pId0 << ", Inserting || spine point "; PrintPoint(outPoint1); std::cout << std::endl;
		spineIds->InsertNextId(pId0);

		if(i < numTrunkPoints - 1)
		{
			spinePointTransform->Translate(stepVector);
		}
	}

	spinePointTransform->RotateWXYZ(-45.0, 0.0, 0.0, 1.0);
	spinePointTransform->Translate(stepVector);

	for(int i = 0; i < numBranchPoints; i++)
	{
		double outPoint1[3];
		spinePointTransform->TransformPoint(rootPoint, outPoint1);
		vtkIdType pId0 = spinePoints->InsertNextPoint(outPoint1);
		std::cout << pId0 << ", Inserting // spine point " << pId0 << ": "; PrintPoint(outPoint1); std::cout << std::endl;
		spineIds->InsertNextId(pId0);

		if(i < numBranchPoints - 1)
		{
			spinePointTransform->Translate(stepVector);
		}
	}
	std::cout << "Total of " << spinePoints->GetNumberOfPoints() << " points in the spine." << std::endl;

	spineLines->InsertNextCell(spineIds);

	vtkSmartPointer<vtkPolyData> polySpine = vtkSmartPointer<vtkPolyData>::New();
	polySpine->SetPoints(spinePoints);
	polySpine->SetLines(spineLines);
	polySpine->GetPointData()->SetScalars(spineRadii);

#if 0
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(polySpine);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(512, 512);
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkInteractorStyleSwitch *iss = vtkInteractorStyleSwitch::SafeDownCast(renderWindowInteractor->GetInteractorStyle());
	iss->SetCurrentStyleToTrackballCamera();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderer->AddActor(actor);

	renderWindow->Render();
	renderWindowInteractor->Start();
	renderWindow->Finalize();
#endif
	// Done with the spine.

	// Insert boundary points. The boundary has four segments.
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();

	int cQuads = 18; // m = 17. Num quads should be even, to make sure m is odd.
	int yQuads = polySpine->GetNumberOfPoints() - 1; // n = 59. Num quads should be even, to make sure n is odd.

	// This is the total number of points to have in the boundary.
	vtkIdType pIds = (cQuads + yQuads) * 2;

	double p0[3]; // Tmp storage.
	double p1[3]; // Tmp storage.

	// TODO: Assuming the polySpline has at least two points. Actually it should have at least 4, as dictated by dbihar.
	double inletDirection[3];
	polySpine->GetPoint(0, p0);
	polySpine->GetPoint(1, p1);
	vtkMath::Subtract(p1, p0, inletDirection);

	// TODO: Assuming the polySpline has at least two points. Actually it should have at least 4, as dictated by dbihar.
	double outletDirection[3];
	polySpine->GetPoint(polySpine->GetNumberOfPoints() - 1, p0);
	polySpine->GetPoint(polySpine->GetNumberOfPoints() - 2, p1);
	vtkMath::Subtract(p1, p0, outletDirection);

	// TODO: This should be provided with the centreline; it should be calculated to match the angle of the cross product
	// of the bifurcation vectors. Perhaps it can be associated with the centreline as vtkInformation object.
	double radiusDirection[3] = {0.0, 0.0, 1.0};

	vtkSmartPointer<vtkTransform> radiusSpinTransform = vtkSmartPointer<vtkTransform>::New();
	vtkSmartPointer<vtkTransform> radiusTranslationTransform = vtkSmartPointer<vtkTransform>::New();

	double outletRotationAngle;
	double outletRotationAxis[3];

	for(vtkIdType pId = 0; pId < pIds; pId++)
	{
		// This point is declared inside the loop to make sure it is (0,0,0) at the start of every iteration.
		double point[3] = {0.0};

		// Terminal arc.
		// Inserting points along the y = y1 boundary segment (in theory).
		double localRadius = polySpine->GetPointData()->GetScalars()->GetComponent(0, 0);
		// TODO: Use vtkMaths::MultiplyScalar(...).

		double radiusVector[3];
		memcpy(radiusVector, radiusDirection, sizeof(double) * 3);
		vtkMath::MultiplyScalar(radiusVector, localRadius);

		if(pId < cQuads)
		{
			// Parametric angle.
			double dA = pId / (double)cQuads;

			// Set up the transform.
			radiusSpinTransform->Identity();
			// Working around the first point (inlet) of the spine.
			radiusSpinTransform->Translate(polySpine->GetPoint(0));
			// Rotate around the vector pointing from the first spine point to the second.
			radiusSpinTransform->RotateWXYZ(vtkMath::DegreesFromRadians(vtkMath::Pi() * dA), inletDirection);

			// Turn the radius vector.
			radiusSpinTransform->TransformPoint(radiusVector, point);
		}
		// Line parallel to the centreline.
		// Inserting points along the x = x2 boundary segment (in theory).
		else if(pId < cQuads + yQuads)
		{
			// Point id along the centreline.
			int lId = pId - cQuads;

			// Get the point from the centreline, update the transform.
			radiusTranslationTransform->Identity();
			radiusTranslationTransform->Translate(polySpine->GetPoint(lId));

			// Translate the radius vector.
			radiusTranslationTransform->TransformPoint(NegateVector(radiusVector), point);
		}
		// Terminal arc.
		// Inserting points along the y = y2 boundary segment.
		else if(pId < cQuads * 2 + yQuads)
		{
			// Local point id.
			int lId = pId - cQuads - yQuads;
			// Parametric angle.
			double dA = lId / (double)cQuads;

			// Set up the transform.
			radiusSpinTransform->Identity();
			// Working around the last point (outlet) of the spine.
			radiusSpinTransform->Translate(polySpine->GetPoint(polySpine->GetNumberOfPoints() - 1));
			// Rotate around the vector pointing from the last spine point to the one before that.
			radiusSpinTransform->RotateWXYZ(vtkMath::DegreesFromRadians(vtkMath::Pi() * dA), outletDirection);

			// Turn the radius vector.
			radiusSpinTransform->TransformPoint(NegateVector(radiusVector), point);
		}
		// Line parallel to the centreline.
		// Inserting points along the x = x1 boundary segment.
		else
		{
			int lId = -1 * (pId - cQuads - yQuads - cQuads) + (polySpine->GetNumberOfPoints() - 1);

			// Get the point from the centreline, update the transform.
			radiusTranslationTransform->Identity();
			radiusTranslationTransform->Translate(polySpine->GetPoint(lId));

			// Translate the radius vector.
			radiusTranslationTransform->TransformPoint(radiusVector, point);
		}
		std::cout << "Ineserting point: "; PrintPoint(point); std::cout << std::endl;
		vtkIdType id = points->InsertNextPoint(point);
		// Sanity check.
		assert(id == pId);
		boundary->GetPointIds()->InsertNextId(pId);
	}
	boundary->GetPointIds()->InsertNextId(0);

	vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
	boundaries->InsertNextCell(boundary);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(points);
	inputPatch->SetLines(boundaries);

	showPolyData(inputPatch, NULL);

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	// Set the bounds of the UV space.
	patchFilter->SetA(0.0);
	patchFilter->SetB(2.0/3.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(vtkMath::Pi());
	// Set the boundary conditions.
	patchFilter->SetMQuads(cQuads);
	patchFilter->SetNQuads(yQuads);
	// Set solution method.
	patchFilter->SetIFlag(2);

	patchFilter->SetInputData(inputPatch);

	// patchFilter->Print(std::cout);

	patchFilter->Update();

	// patchFilter->Print(std::cout);

	vtkPolyData *outputPatch = patchFilter->GetOutput();

	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetDimensions(cQuads + 1, yQuads + 1, 1);
	structuredGrid->SetPoints(outputPatch->GetPoints());

	showPolyData(inputPatch, structuredGrid);

	std::cout << "Exiting " << __FILE__ << std::endl;

	delete [] radiiValues;

	return EXIT_SUCCESS;
}

