#include <cmath>
#include <algorithm>

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>

#include "vtkDbiharPatchFilter.h"

vtkStandardNewMacro(vtkDbiharPatchFilter);

vtkDbiharPatchFilter::vtkDbiharPatchFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->MDim = 0.0;
	this->NDim = 0.0;
	this->IFlag = 0;

	this->Alpha = 0.0;
	this->Beta = 0.0;
	this->Tol = 0.0;
	this->ITCG = 50;

}

int vtkDbiharPatchFilter::RequestData(vtkInformation *vtkNotUsed(request),
		vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	//input->Print(std::cout);
	//output->Print(std::cout);

	// Lines are in the input vtkPolyData.

	// Always want 4 borders.
	assert(input->GetNumberOfLines() == 4);

	input->GetLines()->InitTraversal();

	vtkSmartPointer<vtkIdList> lowBorder = vtkIdList::New();
	input->GetLines()->GetNextCell(lowBorder);

	vtkSmartPointer<vtkIdList> highBorder = vtkIdList::New();
	input->GetLines()->GetNextCell(highBorder);

	vtkSmartPointer<vtkIdList> leftBorder = vtkIdList::New();
	input->GetLines()->GetNextCell(leftBorder);

	vtkSmartPointer<vtkIdList> rightBorder = vtkIdList::New();
	input->GetLines()->GetNextCell(rightBorder);

	std::cout << lowBorder->GetNumberOfIds() << ", " << highBorder->GetNumberOfIds() << ", " << leftBorder->GetNumberOfIds() << ", " << rightBorder->GetNumberOfIds() << std::endl;

	assert(lowBorder != NULL);
	assert(highBorder != NULL);

	assert(leftBorder != NULL);
	assert(rightBorder != NULL);

	// Low and high borders must have the same number of points.
	assert(lowBorder->GetNumberOfIds() == highBorder->GetNumberOfIds());

	// Low and high borders must have the same number of points.
	assert(leftBorder->GetNumberOfIds() == rightBorder->GetNumberOfIds());

	this->MDim = lowBorder->GetNumberOfIds() - 2;
	this->NDim = leftBorder->GetNumberOfIds() - 2;

	this->Print(std::cout);

	vtkPoints *outputPoints = vtkPoints::New();

	std::cout << outputPoints->GetNumberOfPoints() << std::endl;

	double tmpPoint[3] = {0.0, 0.0, 0.0};
	for(int row = 0; row < leftBorder->GetNumberOfIds(); row++)
	{
		for(int col = 0; col < lowBorder->GetNumberOfIds(); col++)
		{
			outputPoints->InsertNextPoint(tmpPoint);
		}
	}

	std::cout << outputPoints->GetNumberOfPoints() << std::endl;

	// Translate the centre of the dataset to zero?

	double a = 1.0;
	double b = 3.0;
	double c = 1.0;
	double d = 4.0;

	// Allocate derivatives. And initialise them to 0, at leas for now.
	double *bda = new double[this->NDim];
	std::fill_n(bda, this->NDim, 0.0);
	double *bdb = new double[this->NDim];
	std::fill_n(bdb, this->NDim, 0.0);
	double *bdc = new double[this->MDim];
	std::fill_n(bdc, this->MDim, 0.0);
	double *bdd = new double[this->MDim];
	std::fill_n(bdd, this->MDim, 0.0);

	// Allocate f;
	double **f = new double*[this->NDim + 2];
	for(int i = 0; i < this->NDim + 2; i++)
	{
		f[i] = new double[this->MDim + 2];
	}

	int idf = this->MDim + 2;

	// From the description of Dbihar source code in Fortran.
	int lw = std::max(3 * this->MDim, 4 * this->NDim) + 4 * this->NDim + 2 * this->MDim +0.5 * pow(this->MDim + 1, 2) + 19;

	// Allocate workspace.
	double *w = new double[lw];

	// Do this once for each dimension, X, Y, Z.
	for(int dim = 0; dim < 3; dim++)
	{
		// Copy points coordinates (per current dimension) into f;
		std::cout << "Inserting low and high boundaries..." << std::endl;
		for(int i = 0; i < this->MDim + 2; i++)
		{
			// Get the coordinate value for this dimension and put it into f along the low boundary.
			double *p1 = input->GetPoint(lowBorder->GetId(i));
			f[0][i] = p1[dim];
			std::cout << 0 << ":" << i << ", [" << p1[0] << "," << p1[1] << "," << p1[2] << "]" << std::endl;
			// Get the coordinate value for this dimension and put it f along the high boundary.
			double *p2 = input->GetPoint(highBorder->GetId(i));
			f[this->NDim + 1][i] = p2[dim];
			std::cout << this->NDim + 1 << ":" << i << ", [" << p2[0] << "," << p2[1] << "," << p2[2] << "]" << std::endl;
		}

		// Copy points coordinates (per current dimension) into f;
		std::cout << "Inserting left and right boundaries..." << std::endl;
		for(int i = 0; i < this->NDim + 2; i++)
		{
			// Get the coordinate value for this dimension and put it into f along the left boundary.
			double *p1 = input->GetPoint(leftBorder->GetId(i));
			f[i][0] = p1[dim];
			std::cout << i << ":" << 0 << ", [" << p1[0] << "," << p1[1] << "," << p1[2] << "]" << std::endl;
			// Get the coordinate value for this dimension and put it f along the right boundary.
			double *p2 = input->GetPoint(rightBorder->GetId(i));
			f[i][this->MDim + 1] = p2[dim];
			std::cout << i << ":" << this->MDim + 1 << ", [" << p2[0] << "," << p2[1] << "," << p2[2] << "]" << std::endl;
		}

		for(int n = 0; n < this->NDim + 2; n++)
		{
			for(int m = 0; m < this->MDim + 2; m++)
			{
				std::cout << f[n][m] << " ";
			}
			std::cout << std::endl;
		}

		dbihar_(&a, &b, &(this->MDim),
				bda, bdb, bdc, bdd,
				&c, &d, &(this->NDim),
				(double **)f, &idf,
				&(this->Alpha), &(this->Beta), &(this->IFlag), &(this->Tol), &(this->ITCG),
				w, &lw);

		std::cout << "Return IFlag: " << this->IFlag << std::endl;

		// Save result.
		vtkIdType pId = 0;
		for(int row = 0; row < leftBorder->GetNumberOfIds(); row++)
		{
			for(int col = 0; col < lowBorder->GetNumberOfIds(); col++, pId++)
			{
				outputPoints->GetPoint(pId, tmpPoint);
				std::cout << pId << ", " << row << ":" << col << ", [" << tmpPoint[0] << "," << tmpPoint[1] << "," << tmpPoint[2] << "]" << std::endl;
				tmpPoint[dim] = f[row][col];
				std::cout << "[" << tmpPoint[0] << "," << tmpPoint[1] << "," << tmpPoint[2] << "]" << std::endl;
				outputPoints->InsertPoint(pId, tmpPoint);
			}
		}

		// All done for this dimension.
	}

	//output->ShallowCopy(input);
	output->SetPoints(outputPoints);
	std::cout << output->GetPoints()->GetNumberOfPoints() << std::endl;

	// Deallocate f;
	for(int i = 0; i < this->NDim + 2; i++)
	{
		delete [] f[i];
	}
	delete [] f;

	// Deallocated derivatives.
	delete [] bda;
	delete [] bdb;
	delete [] bdc;
	delete [] bdd;

	// Deallocate workspace.
	delete [] w;

	// Required to return 1 by VTK API.
	return 1;
}

void vtkDbiharPatchFilter::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);

	os << indent << "MDim: " << this->MDim << "\n";
	os << indent << "NDim: " << this->NDim << "\n";
	os << indent << "IFlag: " << this->IFlag << "\n";
}
