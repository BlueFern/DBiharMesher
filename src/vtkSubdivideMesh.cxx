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

	if (this->Columns <= 0 || this->Rows <= 0)
	{
		vtkErrorMacro("Both Columns and Rows must be positive numbers.");
		exit(EXIT_FAILURE);
	}

	int numberOfQuads = input->GetNumberOfCells();

	vtkSmartPointer<vtkIntArray> branchIdArray = vtkIntArray::SafeDownCast(input->GetCellData()->GetArray(vtkDbiharStatic::BRANCH_ID_ARR_NAME));
	vtkSmartPointer<vtkIntArray> gridCoordsArray = vtkIntArray::SafeDownCast(input->GetCellData()->GetArray(vtkDbiharStatic::GRID_COORDS_ARR_NAME));

	std::map<int, std::pair<int, int> > branchDimsensions;
	std::map<int, std::pair<int, int> >::iterator iter;

	// For each branch retrieve the highest number of quad rows and columns in the grid.
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
				iter->second.first = rowVal;
			if(iter->second.second < colVal)
				iter->second.second = colVal;
		}
	}

	std::map<int, int> labelRowOffsets;

	// We need to remember the maximum number of rows per label.
	// Like a cumulative value.
	std::map<int, int> maxRowsPerLabel;
	for(iter = branchDimsensions.begin(); iter != branchDimsensions.end(); iter++)
	{
		maxRowsPerLabel[iter->first] = 0;
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
		// TODO: What happened to the good old loop?
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

		// TODO: Add a bool flag to this filter which disables the parametric grid coordinates calculation by default.

		// For the processed quad add grid coordinates array.
		// Iterate over the cells of the subdivided quad and on the basis
		// of the number of current quad and the number of quads in each row,
		// calculate grid coordinates.
		int branchId = branchIdArray->GetValue(quadId);
		int numberOfCellsPerQuad = subdividedQuad->GetNumberOfCells();

		vtkSmartPointer<vtkIntArray> gridCoordinatesCellData = vtkSmartPointer<vtkIntArray>::New();
		gridCoordinatesCellData->SetName(vtkDbiharStatic::GRID_COORDS_ARR_NAME);
		gridCoordinatesCellData->SetNumberOfComponents(2);

		for(int cellId = 0; cellId < numberOfCellsPerQuad; cellId++)
		{
			// Figure out the column and row.
			int rowId = cellId / this->Rows;
			int colId = cellId - rowId * this->Rows;

			// Remember the first global row id as offset for this label.
			int globalRowId = (quadId / branchDimsensions[branchId].second) * this->Columns + rowId;
			std::map<int, int>::iterator iter = labelRowOffsets.find(branchId);
			if(iter == labelRowOffsets.end())
				labelRowOffsets[branchId] = globalRowId;

			// Roll back the row id to start from zero for this label.
			gridCoords[0] = globalRowId - labelRowOffsets[branchId];
			gridCoords[1] = (quadId % branchDimsensions[branchId].second) * this->Rows + colId;

			gridCoordinatesCellData->InsertNextTuple(gridCoords);
		}

		subdividedQuad->GetCellData()->AddArray(gridCoordinatesCellData);

		// TODO: Add branch id array for all cells for future use?

		appendPolyData->AddInputData(subdividedQuad);

		// Basic progress reporting.
		if ((quadId + 1) % percentOfTotal == 0)
		{
			// 10% + 1 as we start with stage at 1.
			this->UpdateProgress(static_cast<double>(stage++) / static_cast<double>(11));
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
