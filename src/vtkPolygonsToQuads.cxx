#include <algorithm>

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkCallbackCommand.h>
#include <vtkSmartPointer.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>

#include "vtkPolygonsToQuads.h"
#include "vtkDbiharStatic.h"
#include "vtkSetGet.h"

vtkStandardNewMacro(vtkPolygonsToQuads);

vtkPolygonsToQuads::vtkPolygonsToQuads()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

}

int vtkPolygonsToQuads::RequestData(vtkInformation *vtkNotUsed(request),
		vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();

	vtkSmartPointer<vtkCellData> oldCellData = input->GetCellData();
	int numberOfArrays = oldCellData->GetNumberOfArrays();

	vtkSmartPointer<vtkDoubleArray> oldCellDataArrays[numberOfArrays];
	vtkSmartPointer<vtkDoubleArray> newCellDataArrays[numberOfArrays];

	// Get all previous cell data and initialise new cell data arrays with their
	// appropriate names.
	for (int i = 0; i < numberOfArrays; i++)
	{
		oldCellDataArrays[i] = vtkDoubleArray::SafeDownCast(oldCellData->GetArray(i));
		newCellDataArrays[i] = vtkSmartPointer<vtkDoubleArray>::New();
		newCellDataArrays[i]->SetName(oldCellDataArrays[i]->GetName());
	}

	vtkSmartPointer<vtkIdList> oldCell = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkIdList> newCell = vtkSmartPointer<vtkIdList>::New();
	newCell->SetNumberOfIds(4);

	// Iterate over all cells. Create two quads from all 6-sided polygons and give both
	// new quads the same cell data as the original polygon.
	for (int i = 0; i < input->GetNumberOfCells(); i++)
	{
		if (input->GetCell(i)->GetNumberOfPoints() == 6)
		{
			input->GetCellPoints(i, oldCell);

			newCell->SetId(0, oldCell->GetId(0));
			newCell->SetId(1, oldCell->GetId(1));
			newCell->SetId(2, oldCell->GetId(4));
			newCell->SetId(3, oldCell->GetId(5));
			cells->InsertNextCell(newCell);

			newCell->SetId(0, oldCell->GetId(1));
			newCell->SetId(1, oldCell->GetId(2));
			newCell->SetId(2, oldCell->GetId(3));
			newCell->SetId(3, oldCell->GetId(4));
			cells->InsertNextCell(newCell);

			for (int j = 0; j < numberOfArrays; j++)
			{
				newCellDataArrays[j]->InsertNextTuple(oldCellDataArrays[j]->GetTuple(i));
				newCellDataArrays[j]->InsertNextTuple(oldCellDataArrays[j]->GetTuple(i));
			}
		}

		// Already a quad, copy to new cell array and copy the cell data for the quad.
		else
		{
			cells->InsertNextCell(input->GetCell(i));

			for (int j = 0; j < numberOfArrays; j++)
			{
				newCellDataArrays[j]->InsertNextTuple(oldCellDataArrays[j]->GetTuple(i));
			}
		}
	}

	// Create new polydata with new cells and cell data.
	polydata->SetPoints(input->GetPoints());
	polydata->SetPolys(cells);

	for (int j = 0; j < numberOfArrays; j++)
	{
		polydata->GetCellData()->AddArray(newCellDataArrays[j]);
	}
	output->ShallowCopy(polydata);

	return 1;
}

