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

	this->MDim = 0;
	this->NDim = 0;
	this->MQuads = 0;
	this->NQuads = 0;
	this->IFlag = 0;
	this->OFlag = 0;

	this->Alpha = 0.0;
	this->Beta = 0.0;
	this->Tol = 1e-3;
	this->ITCG = 10;
}

int vtkDbiharPatchFilter::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	// Test the number of boundary points matches the expected number of
	// points according to the passed values of MQuads and NQuads.
	// TODO: Provide a helpful output message if this fails.
	vtkIdType pIds = (this->MQuads + this->NQuads) * 2;
	assert(input->GetNumberOfPoints() == pIds);

	// TODO: Review whether MDim and NDim members can be removed.
	this->MDim = this->MQuads - 1;
	this->NDim = this->NQuads - 1;

	// Prepare output points.
	vtkPoints *outputPoints = vtkPoints::New();
	double tmpPoint[3] = {0.0, 0.0, 0.0};
	for(int n = 0; n < (this->NQuads + 1); n++)
	{
		for(int m = 0; m < (this->MQuads + 1); m++)
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
		for(vtkIdType pId = 0; pId < pIds; pId++)
		{

			double val = input->GetPoint(pId)[dim];
			int fIdx = 0;
			// Inserting from the y = y1 boundary segment.
			if(pId < this->MQuads)
			{
				fIdx = pId;
			}
			// Inserting from the x = x2 boundary segment.
			else if(pId < this->MQuads + this->NQuads)
			{
				int locId = pId - this->MQuads;
				fIdx = locId * (this->MQuads + 1) + this->MQuads;
			}
			// Inserting from the y = y2 boundary segment.
			else if(pId < this->MQuads * 2 + this->NQuads)
			{
				int locId = abs((pId - this->MQuads - this->NQuads) - this->MQuads);
				fIdx = locId + (this->NQuads * (this->MQuads + 1));
			}
			// Inserting from the x = x1 boundary segment.
			else
			{
				int locId = abs((pId - this->MQuads * 2 - this->NQuads) - this->NQuads);
				fIdx = locId * (this->MQuads + 1);
			}
			f[fIdx] = val;
		}

#if PRINT_DEBUG
		int ind = 0;
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
#if 1
		if(dim == 1)
		{
			std::fill_n(bda, this->NDim, 0.0);
			std::fill_n(bdb, this->NDim, 0.0);
			std::fill_n(bdc, this->MDim, -10.0);
			std::fill_n(bdd, this->MDim, 10.0);
		}

		if(dim == 2)
		{
			std::fill_n(bda, this->NDim, -60.0);
			std::fill_n(bdb, this->NDim, -60.0);
			std::fill_n(bdc, this->MDim, 0.0);
			std::fill_n(bdd, this->MDim, 0.0);
		}
#endif

		this->OFlag = this->IFlag;

		dbihar_(&(this->A), &(this->B), &(this->MDim),
				bda, bdb, bdc, bdd,
				&(this->C), &(this->D), &(this->NDim),
				(double *)f, &idf,
				&(this->Alpha), &(this->Beta), &(this->OFlag), &(this->Tol), &(this->ITCG),
				w, &lw);

		std::cout << "Initial IFlag: " << this->IFlag << std::endl;
		std::cout << "Returned OFlag: " << this->OFlag << std::endl;
		std::cout << "Returned Tol: " << this->Tol << std::endl;
		std::cout << "Returned ITCG: " << this->ITCG << std::endl;

#if PRINT_DEBUG
		ind = 0;
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
		for(int row = 0; row < this->NQuads + 1; row++)
		{
			for(int col = 0; col < this->MQuads + 1; col++, pId++)
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
	os << indent << "MQuads: " << this->MQuads << "\n";
	os << indent << "NQuads: " << this->NQuads << "\n";

	os << indent << "IFlag: " << this->IFlag << "\n";
	os << indent << "OFlag: " << this-OFlag << "\n";

	os << indent << "Number of input points: " << vtkPolyData::SafeDownCast(this->GetInput())->GetNumberOfPoints() << "\n";
	os << indent << "Number of output points: " << vtkPolyData::SafeDownCast(this->GetOutput())->GetNumberOfPoints() << "\n";

	//os << indent << "Input:" << "\n";
	//this->GetInput()->PrintSelf(os, indent.GetNextIndent());
	//os << indent << "Output:" << "\n";
	//this->GetOutput()->PrintSelf(os, indent.GetNextIndent());
}
