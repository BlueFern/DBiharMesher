/*
 * Circullar Dbihar Patch Generation.
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

#include "vtkDbiharPatchFilter.h"
#include "showPolyData.h"

#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	// Bulid a circular patch.

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

	vtkSmartPointer<vtkCellArray> boundaries = vtkSmartPointer<vtkCellArray>::New();
	boundaries->InsertNextCell(boundary);

	vtkSmartPointer<vtkPolyData> inputPatch = vtkSmartPointer<vtkPolyData>::New();
	inputPatch->SetPoints(points);
	inputPatch->SetLines(boundaries);

	vtkSmartPointer<vtkXMLPolyDataWriter> inputPatchWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	inputPatchWriter->SetFileName("inputPatch.vtp");
	inputPatchWriter->SetInputData(inputPatch);
	inputPatchWriter->Update();

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

	showPolyData(inputPatch, structuredGrid);

	vtkSmartPointer<vtkStructuredGridGeometryFilter> gridGeometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	gridGeometryFilter->SetInputData(structuredGrid);
	gridGeometryFilter->Update();

	showPolyData(gridGeometryFilter->GetOutput(), NULL);

	vtkSmartPointer<vtkXMLPolyDataWriter> polyDataWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	polyDataWriter->SetInputData(gridGeometryFilter->GetOutput());
	polyDataWriter->SetFileName("endCapQuads.vtp");
	polyDataWriter->Update();

	vtkSmartPointer<vtkTriangleFilter> triangulatorFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangulatorFilter->SetInputData(gridGeometryFilter->GetOutput());
	triangulatorFilter->Update();

	showPolyData(triangulatorFilter->GetOutput(), NULL);

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
