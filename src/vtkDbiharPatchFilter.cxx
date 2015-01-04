/*
 * Program: vtkDbiharPatchFilter.
 */

#include <cmath>
#include <algorithm>

#include <vtkObjectFactory.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>

#include <vtkPolyLine.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCallbackCommand.h>

#include "vtkDbiharPatchFilter.h"

#define PRINT_DEBUG 0

vtkStandardNewMacro(vtkDbiharPatchFilter);

const char *vtkDbiharPatchFilter::DERIV_ARR_NAME = {"derivVectors"};

vtkDbiharPatchFilter::vtkDbiharPatchFilter()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	OutputPoints = vtkSmartPointer<vtkPoints>::New();

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

	vtkSmartPointer<vtkCallbackCommand> progressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	progressCallback->SetCallback(this->ProgressFunction);
	this->AddObserver(vtkCommand::ProgressEvent, progressCallback);
}

int vtkDbiharPatchFilter::RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	// Test the number of boundary points matches the expected number of
	// points according to the passed values of MQuads and NQuads.
	vtkIdType pIds = (this->MQuads + this->NQuads) * 2;
	if(input->GetNumberOfPoints() != pIds)
	{
		vtkErrorMacro("Number of points in the input data (" << input->GetNumberOfPoints() << ") does not match the expected number of points (" << pIds << ").");
		exit(EXIT_FAILURE);
	}

	// TODO: Review whether MDim and NDim members can be removed.
	this->MDim = this->MQuads - 1;
	this->NDim = this->NQuads - 1;

	// Allocate derivatives.
	double *bda = new double[this->NDim];
	//std::fill_n(bda, this->NDim, 0.0);
	double *bdb = new double[this->NDim];
	//std::fill_n(bdb, this->NDim, 0.0);
	double *bdc = new double[this->MDim];
	//std::fill_n(bdc, this->MDim, 0.0);
	double *bdd = new double[this->MDim];
	//std::fill_n(bdd, this->MDim, 0.0);

	vtkDataArray *derivatives = input->GetPointData()->GetVectors(DERIV_ARR_NAME);
	if(derivatives == 0)
	{
		vtkWarningMacro("Boundary derivatives are not set in the input data. Setting derivatives to zero.");
		std::fill_n(bda, this->NDim, 0.0);
		std::fill_n(bdb, this->NDim, 0.0);
		std::fill_n(bdc, this->MDim, 0.0);
		std::fill_n(bdd, this->MDim, 0.0);
	}
	else if(derivatives->GetNumberOfTuples() != pIds)
	{
		vtkErrorMacro("Number of derivative vectors in the input data (" << derivatives->GetNumberOfTuples() << ") does not match the expected number of points (" << pIds << ").");
		exit(EXIT_FAILURE);
	}

	// Allocate f.
	double *f = new double[(this->NDim + 2) * (this->MDim + 2)];
	// Dbihar require this.
	int idf = this->MDim + 2;

	// From the description of Dbihar source code in Fortran.
	int lw;
	if(this->IFlag == 2)
	{
		// This is the lower bound (?) estimate.
		// lw = (int)(std::max(7 * this->NDim, 3 * this->MDim) + 2 * (this->NDim + this->MDim) + 19);
		// This is the upper bound (?) estimate.
		lw = (int)(std::max(7 * this->NDim, 3 * this->MDim) + 2 * (this->NDim + this->MDim) + 19 + 20 * (this->NDim + 3));
	}
	else if(this->IFlag == 4)
	{
		lw = (int)(std::max(3 * this->MDim, 4 * this->NDim) + 4 * this->NDim + 2 * this->MDim +0.5 * pow(this->NDim + 1, 2) + 19);
	}
	else
	{
		// Other values for IFlag not supported.
		vtkErrorMacro("Unsupported value for IFlag: " << this->IFlag << ".");
		exit(EXIT_FAILURE);
	}

	// Allocate workspace.
	double *w = new double[lw];

	// Prepare output points.
	double tmpPoint[3] = {0.0, 0.0, 0.0};
	for(int pId = 0; pId < (this->NQuads + 1) * (this->MQuads + 1); pId++)
	{
		OutputPoints->InsertNextPoint(tmpPoint);
	}

	int numDims = 3;
	// Do this once for each dimension, X, Y, Z.
	for(int dim = 0; dim < numDims; dim++)
	{
		// Reset f.
		std::fill_n(f, (this->NDim + 2) * (this->MDim + 2), 0.0);

		// Copy points coordinates (per current dimension) into f.
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

		// Prepare derivatives arrays.
		if(derivatives != 0)
		{
			// Copy derivative values (per current dimension) into the appropriate derivative array.
			for(vtkIdType pId = 0; pId < pIds; pId++)
			{
				double val = derivatives->GetComponent(pId, dim);

				// Derivative arrays bda, bdb, bdc, bdd don't have space for the
				// corner points, hence the corner point derivative values are skipped.

				// Inserting derivatives from the y = y1 segment, skipping the corner point.
				if(pId < this->MQuads)
				{
					if(pId != 0)
					{
						int locId = pId - 1;
						assert(locId >= 0 && locId < this->MDim);
						bdc[locId] = val;
					}
				}
				// Inserting derivatives from the x = x2 segment, skipping the corner point.
				else if(pId < this->MQuads + this->NQuads)
				{
					if(pId != this->MQuads)
					{
						int locId = pId - this->MQuads - 1;
						assert(locId >= 0 && locId < this->NDim);
						bdb[locId] = val;
					}
				}
				// Inserting derivatives from the y = y2 segment, skipping the corner point.
				else if(pId < this->MQuads * 2 + this->NQuads)
				{
					if(pId != this->MQuads + this->NQuads)
					{
						int locId = abs((pId - this->MQuads - this->NQuads) - this->MQuads) - 1;
						assert(locId >= 0 && locId < this->MDim);
						bdd[locId] = val;
					}
				}
				// Inserting derivatives from the x = x1  segment, skipping the corner point.
				else
				{
					if(pId != this->MQuads * 2 + this->NQuads)
					{
						int locId = abs((pId - this->MQuads * 2 - this->NQuads) - this->NQuads) - 1;
						assert(locId >= 0 && locId < this->NDim);
						bda[locId] = val;
					}
				}
			}
		}

		this->OFlag = this->IFlag;

		dbihar_(&(this->A), &(this->B), &(this->MDim),
				bda, bdb, bdc, bdd,
				&(this->C), &(this->D), &(this->NDim),
				(double *)f, &idf,
				&(this->Alpha), &(this->Beta), &(this->OFlag), &(this->Tol), &(this->ITCG),
				w, &lw);

		CheckError();

		// std::cout << "Initial IFlag: " << this->IFlag << std::endl;
		// std::cout << "Returned OFlag: " << this->OFlag << std::endl;
		// std::cout << "Returned Tol: " << this->Tol << std::endl;
		// std::cout << "Returned ITCG: " << this->ITCG << std::endl;

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

		// Save result to output points.
		vtkIdType pId = 0;
		// For each element of f copy the value into the current dimension of the corresponding output point.
		for(int row = 0; row < this->NQuads + 1; row++)
		{
			for(int col = 0; col < this->MQuads + 1; col++, pId++)
			{
				OutputPoints->GetPoint(pId, tmpPoint);
				tmpPoint[dim] = f[pId];
				OutputPoints->InsertPoint(pId, tmpPoint);
			}
		}

		if(dim < numDims - 1)
		{
			// Basic progress reporting.
			this->UpdateProgress(static_cast<double>(dim + 1)/static_cast<double>(numDims));
		}

		// All done for this dimension.
	}

	output->SetPoints(OutputPoints);

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

void vtkDbiharPatchFilter::ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
	vtkDbiharPatchFilter* filter = static_cast<vtkDbiharPatchFilter *>(caller);
	cout << filter->GetClassName() << " progress: " << std::fixed << std::setprecision(3) << filter->GetProgress() << endl;
}

void vtkDbiharPatchFilter::CheckError()
{
	switch(this->OFlag)
	{
		case 0:
			vtkErrorMacro("Something is rotten in the state of Denmark, because zero return from dbihar is unexpected.");
			exit(EXIT_FAILURE);
			//return;
		case -1:
			vtkErrorMacro("Dbihar: n and/or m is even or less than 3.");
			exit(EXIT_FAILURE);
			//return;
		case -2:
			vtkErrorMacro("Dbihar: a >= b and/or c >= d.");
			exit(EXIT_FAILURE);
			//return;
		case -3:
			vtkErrorMacro("Dbihar: idf < m + 2 or lw is too small.");
			exit(EXIT_FAILURE);
			//return;
		case -4:
			vtkErrorMacro("Dbihar: Linpack failure in cholesky-factorization. This should not occur, check input carefully.");
			exit(EXIT_FAILURE);
			//return;
		case -5:
			vtkErrorMacro("Dbihar: Linpack detected a computationally singular system using the symmetric indefinite factorization.");
			exit(EXIT_FAILURE);
			//return;
		case -6:
			vtkErrorMacro("Dbihar: The conjugate gradient iteration failed to converge in 30 iterations. The probable cause is an indefinite or near singular system. Try using iflag=4. Note that tol returns an estimate of the residual in the current conjugate gradient iteration.");
			exit(EXIT_FAILURE);
			//return;
		default:
			// TODO: Is this normal completion here?
			;
	}
}
