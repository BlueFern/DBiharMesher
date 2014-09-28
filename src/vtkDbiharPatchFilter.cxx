#include <cmath>
#include <algorithm>

#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>
#include <vtkPolyLine.h>
#include <vtkCellArray.h>

#include "vtkDbiharPatchFilter.h"

#define PRINT_DEBUG 0

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
	this->Tol = 1e-3;
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

	// TODO: Implement better error reporting.
	// Low and high borders must have the same number of points.
	assert(lowBorder->GetNumberOfIds() > 5 && lowBorder->GetNumberOfIds() == highBorder->GetNumberOfIds());

	// TODO: Implement better error reporting.
	// Low and high borders must have the same number of points.
	assert(leftBorder->GetNumberOfIds() > 5 && leftBorder->GetNumberOfIds() == rightBorder->GetNumberOfIds());

	this->MDim = lowBorder->GetNumberOfIds() - 2;
	this->NDim = leftBorder->GetNumberOfIds() - 2;

	vtkPoints *outputPoints = vtkPoints::New();

	double tmpPoint[3] = {0.0, 0.0, 0.0};

	// Prepare output points.
	for(int row = 0; row < leftBorder->GetNumberOfIds(); row++)
	{
		for(int col = 0; col < lowBorder->GetNumberOfIds(); col++)
		{
			outputPoints->InsertNextPoint(tmpPoint);
		}
	}

	// TODO: Derivatives are to be set as parameters of this filter.
	// Allocate derivatives. And initialise them to 0, at least for now.
	double *bda = new double[this->NDim];
	std::fill_n(bda, this->NDim, 0.0);
	double *bdb = new double[this->NDim];
	std::fill_n(bdb, this->NDim, 0.0);
	double *bdc = new double[this->MDim];
	std::fill_n(bdc, this->MDim, 0.0);
	double *bdd = new double[this->MDim];
	std::fill_n(bdd, this->MDim, 0.0);

	// Allocate f.
	double *f = new double[(this->NDim + 2) * (this->MDim + 2)];
	// Dbihar require this.
	int idf = this->MDim + 2;

	// From the description of Dbihar source code in Fortran.
	int lw;
	if(this->IFlag == 2)
	{
		lw = (int)(std::max(7 * this->NDim, 3 * this->MDim) + 2 * (this->NDim + this->MDim) + 19);
	}
	else if(this->IFlag == 4)
	{
		lw = (int)(std::max(3 * this->MDim, 4 * this->NDim) + 4 * this->NDim + 2 * this->MDim +0.5 * pow(this->NDim + 1, 2) + 19);
	}
	else
	{
		// TODO: Provide a better error reporting mechanism.
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
		for(int i = 0; i < this->MDim + 2; i++)
		{
			// Get the coordinate value for this dimension and put it into f along the low boundary.
			input->GetPoint(lowBorder->GetId(i), tmpPoint);
			f[i] = tmpPoint[dim];
			// Get the coordinate value for this dimension and put it f along the high boundary.
			input->GetPoint(highBorder->GetId(i), tmpPoint);
			f[(this->MDim + 2) * (this->NDim + 1) + i] = tmpPoint[dim];
		}

		// Copy points coordinates (per current dimension) into f;
		for(int i = 0; i < this->NDim + 2; i++)
		{
			// Get the coordinate value for this dimension and put it into f along the left boundary.
			input->GetPoint(leftBorder->GetId(i), tmpPoint);
			f[i * (this->MDim + 2)] = tmpPoint[dim];
			// Get the coordinate value for this dimension and put it f along the right boundary.
			input->GetPoint(rightBorder->GetId(i), tmpPoint);
			f[i * (this->MDim + 2) + this->MDim + 1] = tmpPoint[dim];
		}

#if PRINT_DEBUG
		int ind = 0;
		if(this->GetDebug())
		{
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
		}
#endif

		if(dim == 3)
		{
			std::fill_n(bda, this->NDim, 0.0);
			std::fill_n(bdb, this->NDim, 0.0);
			std::fill_n(bdc, this->MDim, 0.0);
			std::fill_n(bdd, this->MDim, 0.0);
		}

		this->OFlag = this->IFlag;

		dbihar_(&(this->A), &(this->B), &(this->MDim),
				bda, bdb, bdc, bdd,
				&(this->C), &(this->D), &(this->NDim),
				(double *)f, &idf,
				&(this->Alpha), &(this->Beta), &(this->OFlag), &(this->Tol), &(this->ITCG),
				w, &lw);

		std::cout << "Returned IFlag: " << this->IFlag << std::endl;
		std::cout << "Returned Tol: " << this->Tol << std::endl;
		std::cout << "Returned ITCG: " << this->ITCG << std::endl;

#if PRINT_DEBUG
		ind = 0;
		if(this->GetDebug())
		{
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
		}
#endif

		// Save result.
		vtkIdType pId = 0;
		// For each element of f copy the value into the current dimension of the corresponding output point.
		for(int row = 0; row < leftBorder->GetNumberOfIds(); row++)
		{
			for(int col = 0; col < lowBorder->GetNumberOfIds(); col++, pId++)
			{
				outputPoints->GetPoint(pId, tmpPoint);
				tmpPoint[dim] = f[col * (this->NDim + 2) + row];
				outputPoints->InsertPoint(pId, tmpPoint);
			}
		}
		// All done for this dimension.
	}

	output->SetPoints(outputPoints);

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

	os << indent << "A: " << this->A << "\n";
	os << indent << "B: " << this->B << "\n";
	os << indent << "C: " << this->C << "\n";
	os << indent << "D: " << this->D << "\n";

	os << indent << "Alpha: " << this->Alpha << "\n";
	os << indent << "Beta: " << this->Beta << "\n";
	os << indent << "Tol: " << this->Tol << "\n";
	os << indent << "ITCG: " << this->ITCG << "\n";

	os << indent << "MDim: " << this->MDim << "\n";
	os << indent << "NDim: " << this->NDim << "\n";
	os << indent << "IFlag: " << this->IFlag << "\n";

	os << indent << "Number of input points: " << vtkPolyData::SafeDownCast(this->GetInput())->GetNumberOfPoints() << "\n";
	os << indent << "Number of output points: " << vtkPolyData::SafeDownCast(this->GetOutput())->GetNumberOfPoints() << "\n";

	//os << indent << "Input:" << "\n";
	//this->GetInput()->PrintSelf(os, indent.GetNextIndent());
	//os << indent << "Output:" << "\n";
	//this->GetOutput()->PrintSelf(os, indent.GetNextIndent());
}
