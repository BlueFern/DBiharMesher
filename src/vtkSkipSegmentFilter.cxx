#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkIntArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCallbackCommand.h>
#include <vtkTransformFilter.h>
#include <vtkAppendPoints.h>
#include <vtkPointSet.h>
#include <vtkTriangleFilter.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkTransform.h>
#include <vtkMath.h>

#include "vtkSkipSegmentFilter.h"
#include "vtkDbiharStatic.h"

vtkStandardNewMacro(vtkSkipSegmentFilter);

vtkSkipSegmentFilter::vtkSkipSegmentFilter()
{
	this->NumberOfRadialQuads = 0;
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
	this->Inlet = false;
	this->Outlet = false;
	this->PointId = -1;
	this->SkipSize = -1;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkSkipSegmentFilter::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	// Check input is good.

	if (this->SkipSize < 0)
	{
		vtkErrorMacro("Positive skip size must be set.");
		exit(EXIT_FAILURE);
	}

	if (this->PointId < 0)
	{
		vtkErrorMacro("Point Id must be set.");
		exit(EXIT_FAILURE);
	}

	if ((this->NumberOfRadialQuads & 1) != 0)
	{
		vtkErrorMacro("Number of radial quads must be even.");
		exit(EXIT_FAILURE);
	}

	if (this->Inlet == this->Outlet)
	{
		vtkErrorMacro("Must set either Inlet or Outlet to true.");
		exit(EXIT_FAILURE);
	}

	// Build ring of points around point specified.
	vtkSmartPointer<vtkDoubleArray> radiiArray =
			vtkDoubleArray::SafeDownCast(input->GetPointData()->GetVectors(vtkDbiharStatic::RADII_VECTORS_ARR_NAME));

	// Temporary buffers.
	double r[3];
	double p0[3];
	double p1[3];
	double v0[3];
	double newTranslation[3] = {0.0};
	double point[3] = {0.0};

	vtkSmartPointer<vtkPoints> patchPoints = vtkSmartPointer<vtkPoints>::New();

	vtkSmartPointer<vtkIdList> lastRing = vtkSmartPointer<vtkIdList>::New();

	// Translation comes from the first point.
	input->GetPoint(this->PointId, p0);

	// Axis of rotation comes from the centreline direction at the first point.
	if (this->Inlet)
	{
		input->GetPoint(this->PointId + 1, p1);
	}
	else
	{
		input->GetPoint(this->PointId - 1, p1);
	}

	vtkMath::Subtract(p0, p1, v0);

	// Get the first radius.
	radiiArray->GetTuple(this->PointId, r);


	for (vtkIdType ptId = 0; ptId < (2 * this->NumberOfRadialQuads); ptId++ )
	{
		double parametricCoord = ptId / (double) this->NumberOfRadialQuads;

		// Angle of rotation comes from the parametric coordinate along the arc.
		double angle = vtkMath::Pi() * parametricCoord;

		// Assemble local transform.
		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		transform->Translate(p0);
		transform->RotateWXYZ(vtkMath::DegreesFromRadians(angle), v0);

		// Transform the radius.
		transform->TransformPoint(r, point);

		patchPoints->InsertNextPoint(point);

	}
	patchPoints->InsertNextPoint(patchPoints->GetPoint(0)); // Duplicate first point to close the loop.

	this->UpdateProgress(static_cast<double>(1) / static_cast<double>(3));

	vtkSmartPointer<vtkAppendPoints> appendPoints = vtkSmartPointer<vtkAppendPoints>::New();
	vtkSmartPointer<vtkPolyData> pointSet = vtkSmartPointer<vtkPolyData>::New();

	pointSet->SetPoints(patchPoints);

	appendPoints->AddInputData(pointSet); // Add the initial ring first.

	// Translate this ring SkipSize times following v0 and add the new points to appendPoints.
	for (int i = 0; i < this->SkipSize - 1; i++)
	{
		vtkMath::Add(v0, newTranslation, newTranslation);

		vtkSmartPointer<vtkTransformFilter> transformPoints = vtkSmartPointer<vtkTransformFilter>::New();
		transformPoints->SetInputData(pointSet);

		vtkSmartPointer<vtkTransform> translateTransform = vtkSmartPointer<vtkTransform>::New();

		translateTransform->Translate(newTranslation);

		transformPoints->SetTransform(translateTransform);
		transformPoints->Update();

		appendPoints->AddInputData(transformPoints->GetOutput());

	}
	appendPoints->Update();

	this->UpdateProgress(static_cast<double>(2) / static_cast<double>(3));

	// Transform poly data containing all points into a structured grid so we can use vtkStructuredGridGeometryFilter to
	// easily build a mesh.

	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetPoints(appendPoints->GetOutput()->GetPoints());
	structuredGrid->SetDimensions(2 * this->NumberOfRadialQuads + 1, this->SkipSize, 1);

	vtkSmartPointer<vtkStructuredGridGeometryFilter> structuredGridGeomFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	structuredGridGeomFilter->SetInputData(structuredGrid);

	int extent[6] = {0, 2 * this->NumberOfRadialQuads + 1, 0, this->SkipSize, 0, 0};
	structuredGridGeomFilter->SetExtent(extent);
	structuredGridGeomFilter->Update();

	// Skip segment requires mesh to be made from triangles.

	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputData(structuredGridGeomFilter->GetOutput());
	triangleFilter->Update();

	output->ShallowCopy(triangleFilter->GetOutput());

	// Build a line around the last ring so a cap can be built later.
	int lastRingStart = (2 * this->NumberOfRadialQuads + 1) * (this->SkipSize - 1);
	vtkSmartPointer<vtkCellArray> line = vtkSmartPointer<vtkCellArray>::New();

	for (vtkIdType pt = lastRingStart; pt < lastRingStart + (2 * this->NumberOfRadialQuads + 1); pt++)
	{
		lastRing->InsertNextId(pt);
	}

	line->InsertNextCell(lastRing);
	output->SetLines(line);

	return 1;
}

void vtkSkipSegmentFilter::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkSkipSegmentFilter* filter = static_cast<vtkSkipSegmentFilter *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}

void vtkSkipSegmentFilter::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Skip Size: " << this->SkipSize << "\n";
	os << indent << "Starting Point ID: " << this->PointId << "\n";
	os << indent << "Number of radial quads: " << this->NumberOfRadialQuads << "\n";
	os << indent << "Inlet: " << this->Inlet << "\n";
	os << indent << "Outlet: " << this->Outlet << "\n";

}
