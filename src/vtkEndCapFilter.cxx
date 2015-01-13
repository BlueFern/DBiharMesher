#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>

#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridGeometryFilter.h>
#include <vtkTriangleFilter.h>

#include "vtkDbiharPatchFilter.h"
#include "vtkEndCapFilter.h"

vtkStandardNewMacro(vtkEndCapFilter);

vtkEndCapFilter::vtkEndCapFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);
}

int vtkEndCapFilter::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	vtkSmartPointer<vtkCellArray> lines = input->GetLines();
	// TODO: Error checking: make sure there is only one line.

	vtkSmartPointer<vtkIdList> line = vtkSmartPointer<vtkIdList>::New();

	lines->GetCell(0, line);

	vtkSmartPointer<vtkPoints> dbhInputPoints = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkPolyLine> dbhInputBoundary = vtkSmartPointer<vtkPolyLine>::New();

	// Extract the boundary points from the input.
	for(int pId = 0; pId < line->GetNumberOfIds() - 1; pId++)
	{
		dbhInputBoundary->GetPointIds()->InsertNextId(dbhInputPoints->InsertNextPoint(input->GetPoint(line->GetId(pId))));
	}
	dbhInputBoundary->GetPointIds()->InsertNextId(0);
	vtkSmartPointer<vtkCellArray> dbhInputBoundaries = vtkSmartPointer<vtkCellArray>::New();
	dbhInputBoundaries->InsertNextCell(dbhInputBoundary);

	// Create input patch for the Dbihar patch filter.
	vtkSmartPointer<vtkPolyData> dbhInputPatch = vtkSmartPointer<vtkPolyData>::New();
	dbhInputPatch->SetPoints(dbhInputPoints);
	dbhInputPatch->SetLines(dbhInputBoundaries);

	dbhInputPatch->Print(std::cout);

	vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();
	// Set the bounds of the UV space.
	patchFilter->SetA(0.0);
	patchFilter->SetB(1.0);
	patchFilter->SetC(0.0);
	patchFilter->SetD(1.0);
	// Set the boundary conditions.
	patchFilter->SetMQuads(dbhInputPatch->GetNumberOfPoints() / 4);
	patchFilter->SetNQuads(dbhInputPatch->GetNumberOfPoints() / 4);
	// Set solution method.
	patchFilter->SetIFlag(2);
	// Set input.
	patchFilter->SetInputData(dbhInputPatch);
	patchFilter->Update();

	// Create a structured grid.
	vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
	grid->SetDimensions(dbhInputPatch->GetNumberOfPoints() / 4 + 1, dbhInputPatch->GetNumberOfPoints() / 4 + 1, 1);
	grid->SetPoints(patchFilter->GetOutput()->GetPoints());

	// Convert the grid to a mesh.
	vtkSmartPointer<vtkStructuredGridGeometryFilter> gridGeometryFilter = vtkSmartPointer<vtkStructuredGridGeometryFilter>::New();
	gridGeometryFilter->SetInputData(grid);

	// Convert the quad mesh to a tri mesh.
	vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter->SetInputConnection(gridGeometryFilter->GetOutputPort());
	triangleFilter->Update();

	output->ShallowCopy(triangleFilter->GetOutput());

	// Required to return 1 by VTK API.
	return 1;
}

void vtkEndCapFilter::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);

	//os << indent << "Input:" << "\n";
	//this->GetInput()->PrintSelf(os, indent.GetNextIndent());
	//os << indent << "Output:" << "\n";
	//this->GetOutput()->PrintSelf(os, indent.GetNextIndent());
}

