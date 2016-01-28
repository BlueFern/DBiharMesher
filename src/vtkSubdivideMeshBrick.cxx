#include <map>

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkInformationVector.h>
#include <vtkIdList.h>
#include <vtkPoints.h>
#include <vtkCellData.h>
#include <vtkCallbackCommand.h>
#include <vtkAppendPolyData.h>
#include <vtkPointSet.h>
#include <vtkTriangleFilter.h>

#include "vtkDbiharStatic.h"
#include "vtkSubdivideMeshBrick.h"
#include "vtkJoinSmcBrickMesh.h"
#include "vtkJoinEcBrickMesh.h"
#include "vtkSubdivideQuadBrickEc.h"
#include "vtkSubdivideQuadBrickSmc.h"

vtkStandardNewMacro(vtkSubdivideMeshBrick);

vtkSubdivideMeshBrick::vtkSubdivideMeshBrick()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Columns = 0;
	this->Rows = 0;
	this->Flat = true;
	this->AxialQuads = 0;
	this->CircQuads = 0;
	this->CellType = -1;
	this->Branches = 3;

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkSubdivideMeshBrick::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);

	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	if (this->Columns <= 0 || this->Rows <= 0)
	{
		vtkErrorMacro("Both Columns and Rows must be positive numbers.");
		exit(EXIT_FAILURE);
	}
	if (this->CellType == -1)
	{
		vtkErrorMacro("Must specify an SMC or an EC.");
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

#if 0
	// We need to remember the maximum number of rows per label.
	// Like a cumulative value.
	std::map<int, int> maxRowsPerLabel;
	for(iter = branchDimsensions.begin(); iter != branchDimsensions.end(); iter++)
	{
		maxRowsPerLabel[iter->first] = 0;
	}
#endif

	// Every 10%, for progress function.
	int percentOfTotal  = numberOfQuads / 10;
	int stage = 1;

	vtkSmartPointer<vtkAppendPolyData> appendPolyData = vtkSmartPointer<vtkAppendPolyData>::New();
	vtkSmartPointer<vtkIdList> pointsList = vtkSmartPointer<vtkIdList>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	pointsList->SetNumberOfIds(4);
	points->SetNumberOfPoints(4);

	int quadsInBranch = this->CircQuads * this->AxialQuads;

	double gridCoords[2] = {0, 0};
	for (int quadId = 0; quadId < numberOfQuads; quadId++)
	{
		input->GetCellPoints(quadId, pointsList);
		vtkSmartPointer<vtkPolyData> pointsData = vtkSmartPointer<vtkPolyData>::New();
		vtkSmartPointer<vtkPolyData> subdividedQuad = vtkSmartPointer<vtkPolyData>::New();

		points->SetPoint(0,input->GetPoint(pointsList->GetId(0)));
		points->SetPoint(1,input->GetPoint(pointsList->GetId(1)));
		points->SetPoint(2,input->GetPoint(pointsList->GetId(2)));
		points->SetPoint(3,input->GetPoint(pointsList->GetId(3)));


		pointsData->SetPoints(points);


		if (this->CellType == vtkDbiharStatic::EC)
		{
			vtkSmartPointer<vtkSubdivideQuadBrickEc> subdivideQuadEc = vtkSmartPointer<vtkSubdivideQuadBrickEc>::New();

			if (quadId < quadsInBranch)
			{
				subdivideQuadEc->SetRotated(0);
			}
			else
			{
				subdivideQuadEc->SetRotated(1);
			}
			subdivideQuadEc->SetInputData(pointsData);
			subdivideQuadEc->SetColumns(this->Columns);
			subdivideQuadEc->SetRows(this->Rows);
			subdivideQuadEc->Update();
			subdividedQuad = subdivideQuadEc->GetOutput();
		}


		if (this->CellType == vtkDbiharStatic::SMC)
		{
			vtkSmartPointer<vtkSubdivideQuadBrickSmc> subdivideQuadSmc = vtkSmartPointer<vtkSubdivideQuadBrickSmc>::New();
			if (quadId % this->CircQuads >= this->CircQuads / 2)
			{
				subdivideQuadSmc->SetRotated(1);
			}

			else
			{
				subdivideQuadSmc->SetRotated(0);
			}
			subdivideQuadSmc->SetInputData(pointsData);
			subdivideQuadSmc->SetColumns(this->Columns);
			subdivideQuadSmc->SetRows(this->Rows);
			subdivideQuadSmc->Update();
			subdividedQuad = subdivideQuadSmc->GetOutput();
		}

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
			int rowId = cellId / this->Columns;
			int colId = cellId - rowId * this->Columns;

			// Remember the first global row id as offset for this label.
			int globalRowId = (quadId / branchDimsensions[branchId].second) * this->Rows + rowId;
			std::map<int, int>::iterator iter = labelRowOffsets.find(branchId);
			if(iter == labelRowOffsets.end())
				labelRowOffsets[branchId] = globalRowId;

			// Roll back the row id to start from zero for this label.
			gridCoords[0] = globalRowId - labelRowOffsets[branchId];
			gridCoords[1] = (quadId % branchDimsensions[branchId].second) * this->Columns + colId;

			gridCoordinatesCellData->InsertNextTuple(gridCoords);
		}

		subdividedQuad->GetCellData()->AddArray(gridCoordinatesCellData);


		appendPolyData->AddInputData(subdividedQuad);

		// Basic progress reporting.
		if ((quadId + 1) % percentOfTotal == 0)
		{
			// 10% + 1 as we start with stage at 1.
			this->UpdateProgress(static_cast<double>(stage++) / static_cast<double>(11));
		}
	}
	appendPolyData->Update();

	vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
	polyData = appendPolyData->GetOutput();

	if (this->CellType == vtkDbiharStatic::EC)
	{
		vtkSmartPointer<vtkJoinEcBrickMesh> joinEcBrickMesh = vtkSmartPointer<vtkJoinEcBrickMesh>::New();
		joinEcBrickMesh->SetInputData(polyData);
		joinEcBrickMesh->SetRows(this->Rows);
		joinEcBrickMesh->SetColumns(this->Columns);
		joinEcBrickMesh->SetAxialQuads(this->AxialQuads);
		joinEcBrickMesh->SetCircQuads(this->CircQuads);
		joinEcBrickMesh->SetBranches(this->Branches);
		joinEcBrickMesh->Update();

		output->ShallowCopy(joinEcBrickMesh->GetOutput());
	}

	else if (this->CellType == vtkDbiharStatic::SMC)
	{
		vtkSmartPointer<vtkJoinSmcBrickMesh> joinSmcBrickMesh = vtkSmartPointer<vtkJoinSmcBrickMesh>::New();
		joinSmcBrickMesh->SetInputData(polyData);
		joinSmcBrickMesh->SetRows(this->Rows);
		joinSmcBrickMesh->SetColumns(this->Columns);
		joinSmcBrickMesh->SetAxialQuads(this->AxialQuads);
		joinSmcBrickMesh->SetCircQuads(this->CircQuads);
		joinSmcBrickMesh->SetFlat(this->Flat);
		joinSmcBrickMesh->SetBranches(this->Branches);
		joinSmcBrickMesh->Update();

		output->ShallowCopy(joinSmcBrickMesh->GetOutput());
	}

	return 1;
}

void vtkSubdivideMeshBrick::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Columns: " << this->Columns << "\n";
	os << indent << "Rows: " << this->Rows << "\n";
}

void vtkSubdivideMeshBrick::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkSubdivideMeshBrick* filter = static_cast<vtkSubdivideMeshBrick *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}
