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

	this->A = 0.0;
	this->B = 0.0;
	this->C = 0.0;
	this->D = 0.0;

	this->MDim = 0.0;
	this->NDim = 0.0;
	this->IFlag = 0;

	this->Alpha = 0.0;
	this->Beta = 0.0;
	this->Tol = 1e-7;
	this->ITCG = 0;

}

int vtkDbiharPatchFilter::RequestData(vtkInformation *vtkNotUsed(request),
		vtkInformationVector **inputVector,
        vtkInformationVector *outputVector)
{
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	// Always expect 4 borders/lines.
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

	// std::cout << lowBorder->GetNumberOfIds() << ", " << highBorder->GetNumberOfIds() << ", " << leftBorder->GetNumberOfIds() << ", " << rightBorder->GetNumberOfIds() << std::endl;

	// Low and high borders must have the same number of points.
	assert(lowBorder->GetNumberOfIds() > 5 && lowBorder->GetNumberOfIds() == highBorder->GetNumberOfIds());

	// Low and high borders must have the same number of points.
	assert(leftBorder->GetNumberOfIds() > 5 && leftBorder->GetNumberOfIds() == rightBorder->GetNumberOfIds());

	this->MDim = lowBorder->GetNumberOfIds() - 2;
	this->NDim = leftBorder->GetNumberOfIds() - 2;

	vtkPoints *outputPoints = vtkPoints::New();

	// std::cout << outputPoints->GetNumberOfPoints() << std::endl;

	double tmpPoint[3] = {0.0, 0.0, 0.0};

	// Prepare output points.
	for(int row = 0; row < leftBorder->GetNumberOfIds(); row++)
	{
		for(int col = 0; col < lowBorder->GetNumberOfIds(); col++)
		{
			outputPoints->InsertNextPoint(tmpPoint);
		}
	}

	// std::cout << outputPoints->GetNumberOfPoints() << std::endl;

	// Translate the centre of the dataset to zero?


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
	double *f = new double[(this->NDim + 2) * (this->MDim + 2)];
	//std::fill_n(f, (this->NDim + 2) * (this->MDim + 2), 0.0);

	int idf = this->MDim + 2;

	// From the description of Dbihar source code in Fortran.
	int lw;
	if(this->IFlag == 2)
	{
		// p.lw = (int) (fmax(7 * p.n, 3 * p.m)) + 2 * (p.n + p.m) + 19;
		lw = (int)(std::max(7 * this->NDim, 3 * this->MDim) + 2 * (this->NDim + this->MDim) + 19);
	}
	else if(this->IFlag == 4)
	{
		// p.lw = (int) (fmax(3 * p.m, 4 * p.n)) + 4 * p.n + 2 * p.m + 0.5 * ((p.n + 1) * (p.n + 1)) + 19;
		lw = (int)(std::max(3 * this->MDim, 4 * this->NDim) + 4 * this->NDim + 2 * this->MDim +0.5 * pow(this->NDim + 1, 2) + 19);
	}
	else
	{
		// Other values for IFlag not supported.
		std::cerr << "Unsupported value for IFlag: " << this->IFlag << std::endl;
		exit(EXIT_FAILURE);
	}

	// Allocate workspace.
	double *w = new double[lw];

	// Do this once for each dimension, X, Y, Z.
	for(int dim = 0; dim < 3; dim++)
	{
		// Reset f.
		std::fill_n(f, (this->NDim + 2) * (this->MDim + 2), 0.0);

		// Copy points coordinates (per current dimension) into f;
		// std::cout << "Inserting low and high boundaries..." << std::endl;
		for(int i = 0; i < this->MDim + 2; i++)
		{
			// Get the coordinate value for this dimension and put it into f along the low boundary.
			double *p1 = input->GetPoint(lowBorder->GetId(i));
			// f[i * (this->NDim + 2)] = p1[dim];
			f[i] = p1[dim];
			// Get the coordinate value for this dimension and put it f along the high boundary.
			double *p2 = input->GetPoint(highBorder->GetId(i));
			// f[(i + 1) * (this->NDim + 2) - 1] = p2[dim];
			f[(this->MDim + 2) * (this->NDim + 1) + i] = p2[dim];
		}

		// Copy points coordinates (per current dimension) into f;
		// std::cout << "Inserting left and right boundaries..." << std::endl;
		for(int i = 0; i < this->NDim + 2; i++)
		{
			// Get the coordinate value for this dimension and put it into f along the left boundary.
			double *p1 = input->GetPoint(leftBorder->GetId(i));
			// f[i] = p1[dim];
			f[i * (this->MDim + 2)] = p1[dim];
			// Get the coordinate value for this dimension and put it f along the right boundary.
			double *p2 = input->GetPoint(rightBorder->GetId(i));
			// f[(this->NDim + 2) * (this->MDim + 1) + i] = p2[dim];
			f[i * (this->MDim + 2) + this->MDim + 1] = p2[dim];
		}

		int ind = 0;
		std::cout << "f" << std::endl;
		for(int n = 0; n < this->NDim + 2; n++)
		{
			for(int m = 0; m < this->MDim + 2; m++, ind++)
			{
				std::cout << std::setw(6) << std::setprecision(2) << std::fixed << f[ind] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		dbihar_(&(this->A), &(this->B), &(this->MDim),
				bda, bdb, bdc, bdd,
				&(this->C), &(this->D), &(this->NDim),
				(double *)f, &idf,
				&(this->Alpha), &(this->Beta), &(this->IFlag), &(this->Tol), &(this->ITCG),
				w, &lw);

		std::cout << "Returned IFlag: " << this->IFlag << std::endl;
		std::cout << "Returned Tol: " << this->Tol << std::endl;
		std::cout << "Returned ITCG: " << this->ITCG << std::endl;

		ind = 0;
		std::cout << "f'" << std::endl;
		for(int n = 0; n < this->NDim + 2; n++)
		{
			for(int m = 0; m < this->MDim + 2; m++, ind++)
			{
				std::cout << std::setw(6) << std::setprecision(2) << std::fixed << f[ind] << " ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;

		// Save result.
		vtkIdType pId = 0;
		for(int row = 0; row < leftBorder->GetNumberOfIds(); row++)
		{
			for(int col = 0; col < lowBorder->GetNumberOfIds(); col++, pId++)
			{
				outputPoints->GetPoint(pId, tmpPoint);
				//std::cout << pId << ", " << row << ":" << col << ", [" << tmpPoint[0] << "," << tmpPoint[1] << "," << tmpPoint[2] << "]" << std::endl;
				tmpPoint[dim] = f[col * (this->NDim + 2) + row];
				//std::cout << col * (this->NDim + 2) + row << ", [" << tmpPoint[0] << "," << tmpPoint[1] << "," << tmpPoint[2] << "]" << std::endl;
				outputPoints->InsertPoint(pId, tmpPoint);
			}
		}
		//std::cout << std::endl;

		// All done for this dimension.
	}

	output->SetPoints(outputPoints);
	std::cout << output->GetPoints()->GetNumberOfPoints() << std::endl;

	// Deallocate f;
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
