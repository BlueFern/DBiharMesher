#include <map>

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkInformationVector.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkCallbackCommand.h>
#include <vtkAppendPolyData.h>

#include "vtkDbiharStatic.h"
#include "vtkSubdivideMesh.h"
#include "vtkSubdivideQuadFilter.h"

vtkStandardNewMacro(vtkSubdivideMesh);

vtkSubdivideMesh::vtkSubdivideMesh()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Columns = 0;
	this->Rows = 0;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkSubdivideMesh::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	if (this->Columns == 0.0 || this->Rows == 0.0)
	{
		vtkErrorMacro("Must set both Columns and Rows to positive numbers.");
		exit(EXIT_FAILURE);
	}

	int numberOfQuads = input->GetNumberOfCells();

	vtkSmartPointer<vtkIntArray> branchIdArray = vtkIntArray::SafeDownCast(input->GetCellData()->GetArray(vtkDbiharStatic::BRANCH_ID_ARR_NAME));
	branchIdArray->Print(std::cout);
	vtkSmartPointer<vtkIntArray> gridCoordsArray = vtkIntArray::SafeDownCast(input->GetCellData()->GetArray(vtkDbiharStatic::GRID_COORDS_ARR_NAME));
	gridCoordsArray->Print(std::cout);

	std::map<int, std::pair<int, int> > branchDimsensions;
	std::map<int, std::pair<int, int> >::iterator iter;

	// For each branch retrieve the number of rows and columns in the grid.
	for(int quadId = 0; quadId < numberOfQuads; quadId++)
	{
		int branchId = branchIdArray->GetValue(quadId);
		int rowVal = gridCoordsArray->GetValue(quadId * 2) + 1;
		int colVal = gridCoordsArray->GetValue(quadId * 2 + 1) + 1;

		iter = branchDimsensions.find(branchId);
		if(iter == branchDimsensions.end())
		{
			branchDimsensions[branchId] = std::pair<int, int>(rowVal, colVal);
		}
		else
		{
			if(iter->second.first < rowVal)
			{
				iter->second.first = rowVal;
			}
			if(iter->second.second < colVal)
			{
				iter->second.second = colVal;
			}
		}
	}

	for(std::map<int, std::pair<int, int> >::iterator iter = branchDimsensions.begin(); iter != branchDimsensions.end(); iter++)
	{
		std::cout << iter->first << ", " << iter->second.first << ", " << iter->second.second << std::endl;
	}

	// Every 10%, for progress function.
	int percentOfTotal  = numberOfQuads / 10;
	int stage = 1;

	vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
	vtkSmartPointer<vtkIdList> pointsList = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	pointsList->SetNumberOfIds(4);
	points->SetNumberOfPoints(4);

	double gridCoords[2] = {0, 0};
	for (int quadId = 0; quadId < numberOfQuads; quadId++)
	{
		input->GetCellPoints(quadId, pointsList);
		// TODO: Use GetCell method to get the whole cell instead of getting the points.
		// Like this: input->GetCell(vtkIdType, vtkGenericCell);
		points->SetPoint(0,input->GetPoint(pointsList->GetId(0)));
		points->SetPoint(1,input->GetPoint(pointsList->GetId(1)));
		points->SetPoint(2,input->GetPoint(pointsList->GetId(2)));
		points->SetPoint(3,input->GetPoint(pointsList->GetId(3)));
		vtkSmartPointer<vtkPolyData> pointsData = vtkSmartPointer<vtkPolyData>::New();
		pointsData->SetPoints(points);

		// TODO: Change vtkSubdivideQuadFilter to accept vtkQuad instead of vtkPolyData?
		vtkSmartPointer<vtkSubdivideQuadFilter> subdivideQuadFilter = vtkSmartPointer<vtkSubdivideQuadFilter>::New();
		subdivideQuadFilter->SetInputData(pointsData);
		subdivideQuadFilter->SetColumns(this->Columns);
		subdivideQuadFilter->SetRows(this->Rows);
		subdivideQuadFilter->Update();

		vtkSmartPointer<vtkPolyData> subdividedQuad = subdivideQuadFilter->GetOutput();

		// For the processed quad add grid coordinates array.
		// Iterate over the cells of the subdivided quad and
		// on the basis of the number of current quad and the number of quads in each row,
		// calculate grid coordinates.
		int labelVal = branchIdArray->GetValue(quadId);
		int numberOfCells = subdividedQuad->GetNumberOfCells();

		vtkSmartPointer<vtkIntArray> gridCoordinatesCellData = vtkSmartPointer<vtkIntArray>::New();
		gridCoordinatesCellData->SetName(vtkDbiharStatic::GRID_COORDS_ARR_NAME);
		gridCoordinatesCellData->SetNumberOfComponents(2);

		for(int cellId = 0; cellId < numberOfCells; cellId++)
		{
			// Figure out the column and row.
			int rowId = cellId / this->Rows;
			int colId = cellId - rowId * this->Rows;

			gridCoords[0] = (quadId / branchDimsensions[labelVal].second) * this->Columns + rowId;
			gridCoords[1] = (quadId % branchDimsensions[labelVal].second) * this->Rows + colId;

			gridCoordinatesCellData->InsertNextTuple(gridCoords);
		}

		subdividedQuad->GetCellData()->AddArray(gridCoordinatesCellData);

		// TODO: Add branch id array for all cells.

		appendPolyData->AddInputData(subdividedQuad);

		// Basic progress reporting.
		if ((quadId + 1) % percentOfTotal == 0)
		{
			this->UpdateProgress(static_cast<double>(stage++) / static_cast<double>(11)); // 10% + 1 as we start with stage at 1.
		}
	}
	appendPolyData->Update();

	output->ShallowCopy(appendPolyData->GetOutput());
	return 1;
}

void vtkSubdivideMesh::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Columns: " << this->Columns << "\n";
	os << indent << "Rows: " << this->Rows << "\n";
}

void vtkSubdivideMesh::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkSubdivideMesh* filter = static_cast<vtkSubdivideMesh *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}
