/*
 * Rectangular Dbihar Patch Generation.
 */

#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSTLWriter.h>

#include "vtkDbiharStatic.h"
#include "vtkDbiharPatchFilter.h"
#include "vtkDbiharStatic.h"
#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

#if 0

	// This is a rectangular patch.

	// Create the sides of a patch as 3D points.
	double x1 = -10.0;
	double x2 = 10.0;
	double y1 = -15.0;
	double y2 = 15.0;
	double z = 0.0;

	int xQuads = 20; // m = 19. Num quads should be even, to make sure m is odd.
	int yQuads = 30; // n = 29. Num quads should be even, to make sure n is odd.

	vtkIdType pIds = (xQuads + yQuads) * 2;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();

	// Insert boundary points. The boundary has four segments.
	// The coordinates of the points are calculated specific to the current boundary segment.
	double point[3] = {0.0, 0.0, 0.0};
	for(vtkIdType pId = 0; pId < pIds; pId++)
	{
		// Inserting along the y = y1 boundary segment.
		if(pId < xQuads)
		{
			point[0] = x1 + pId * ((x2 - x1) / xQuads);
			point[1] = y1;
		}
		// Inserting along the x = x2 boundary segment.
		else if(pId < xQuads + yQuads)
		{
			point[0] = x2;
			point[1] = y1 + (pId - xQuads) * ((y2 - y1) / yQuads);
		}
		// Inserting along the y = y2 boundary segment.
		else if(pId < xQuads * 2 + yQuads)
		{
			point[0] = x2 - (pId - xQuads - yQuads) * ((x2 - x1) / xQuads);
			point[1] = y2;
		}
		// Inserting along the x = x1 boundary segment.
		else
		{
			point[0] = x1;
			point[1] = y2 - (pId - xQuads - yQuads - xQuads ) * ((y2 - y1) / yQuads);
		}

		// Z coordinate does not change in this case.
		point[2] = z;

		vtkIdType id = points->InsertNextPoint(point);
		// Sanity check.
		assert(id == pId);
		boundary->GetPointIds()->InsertNextId(pId);
	}
	boundary->GetPointIds()->InsertNextId(0);

#else

	// This is a circular patch.

	int xQuads = 20; // m = 19. Num quads should be even, to make sure m is odd.
	int yQuads = 20; // n = 19. Num quads should be even, to make sure n is odd.
	int numPoints = (xQuads + yQuads) * 2;
	int numPointsPerSide = numPoints / 4;

	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(numPoints);

	vtkSmartPointer<vtkPolyLine> boundary = vtkSmartPointer<vtkPolyLine>::New();
	boundary->GetPointIds()->SetNumberOfIds(numPoints + 1);

	// Insert boundary points. The boundary has four segments.
	double angleInc = (vtkMath::Pi() / 2.0) / numPointsPerSide;
	double radius = 20.0;
	for(vtkIdType pId_1 = 0, pId_2 = numPointsPerSide * 2 - 1, pId_3 = numPointsPerSide * 2, pId_4 = numPoints - 1;
			pId_1 < numPointsPerSide;
			pId_1++, --pId_2, pId_3++, --pId_4)
	{
		double angle = angleInc * pId_1;

		std::cout << vtkMath::DegreesFromRadians(angle) << std::endl;

		double x = cos(angle) * radius;
		double y = sin(angle) * radius;

		double point[3] = {-x, y, 0.0};

		points->InsertPoint(pId_1, point);
		boundary->GetPointIds()->InsertId(pId_1, pId_1);

		point[0] *= -1;
		points->InsertPoint(pId_2, point);
		boundary->GetPointIds()->InsertId(pId_2, pId_2);

		point[1] *= -1;
		points->InsertPoint(pId_3, point);
		boundary->GetPointIds()->InsertId(pId_3, pId_3);

		point[0] *= -1;
		points->InsertPoint(pId_4, point);
		boundary->GetPointIds()->InsertId(pId_4, pId_4);
	}
	boundary->GetPointIds()->InsertId(numPoints, 0);

	for(vtkIdType pId = 0; pId < boundary->GetPointIds()->GetNumberOfIds(); pId++)
	{
		std::cout << boundary->GetPointIds()->GetId(pId) << " ";
	}
	std::cout << std::endl;
#endif

	vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
	boundaries->InsertNextCell(boundary);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(points);
	inputPatch->SetLines(boundaries);

	// showPolyData(inputPatch, NULL);

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

	// Set the bounds of the UV space.
	patchFilter->SetA(0.0);
	//patchFilter->SetB(2.0/3.0);
	patchFilter->SetB(1.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(1.0);
	// Set the boundary conditions.
	patchFilter->SetMQuads(xQuads);
	patchFilter->SetNQuads(yQuads);
	// Set solution method.
	patchFilter->SetIFlag(2);

	patchFilter->SetInputData(inputPatch);

	// patchFilter->Print(std::cout);

	patchFilter->Update();

	//  patchFilter->Print(std::cout);

	vtkPolyData *outputPatch = patchFilter->GetOutput();

	vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
	structuredGrid->SetDimensions(xQuads + 1, yQuads + 1, 1);
	structuredGrid->SetPoints(outputPatch->GetPoints());

	vtkDbiharStatic::ShowPolyDataWithGrid(inputPatch, structuredGrid);

	vtkSmartPointer<vtkStructuredGridGeometryFilter> gridGeometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	gridGeometryFilter->SetInputData(structuredGrid);
	gridGeometryFilter->Update();

	vtkDbiharStatic::ShowPolyData(gridGeometryFilter->GetOutput());

	vtkSmartPointer<vtkXMLPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	polyDataWriter->SetInputData(gridGeometryFilter->GetOutput());
	polyDataWriter->SetFileName("endCapQuads.vtp");
	polyDataWriter->Update();

	vtkSmartPointer<vtkTriangleFilter> triangulatorFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangulatorFilter->SetInputData(gridGeometryFilter->GetOutput());
	triangulatorFilter->Update();

	vtkDbiharStatic::ShowPolyData(triangulatorFilter->GetOutput());

	polyDataWriter->SetInputData(triangulatorFilter->GetOutput());
	polyDataWriter->SetFileName("endCapTriangles.vtp");
	polyDataWriter->Update();

	vtkSmartPointer<vtkSTLWriter> stlWriter = vtkSmartPointer<vtkSTLWriter>::New();
	stlWriter->SetInputData(triangulatorFilter->GetOutput());
	stlWriter->SetFileName("endCapTriangles.stl");
	stlWriter->Update();

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
