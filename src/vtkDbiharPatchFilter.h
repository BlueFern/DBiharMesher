#ifndef __vtkDbiharPatchFilter_h
#define __vtkDbiharPatchFilter_h

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

extern "C" {
void dbihar_(double* a, double* b, int* m,
		double bda[], double bdb[], double bdc[], double bdd[],
		double* c, double* d, int* n, double *f, int* idf,
		double* alpha, double* beta, int* iflag, double* tol, int* itcg,
		double w[], int* lw);
}

class vtkDbiharPatchFilter : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkDbiharPatchFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkDbiharPatchFilter *New();

	vtkSetMacro(A,double);
	vtkSetMacro(B,double);
	vtkSetMacro(C,double);
	vtkSetMacro(D,double);
	//vtkSetMacro(MDim,int);
	//vtkSetMacro(NDim,int);
	vtkSetMacro(Tol,double);
	vtkSetMacro(IFlag,int);
	vtkSetMacro(Alpha,double);
	vtkSetMacro(Beta,double);

	enum ErrorIds
		: signed {
			NoError = 0,
			InvalidMN = -1, ///< n and/or m was even or less than 3.
			InvalidABCDE = -2, ///< a.ge.b and/or c.ge.d .
			InvalidIdf = -3, ///< idf.lt.m+2 or lw is too small.
			LinpackFail = -4, ///< linpack failure in cholesky-factorization.
							  ///< this should not occur,check input carefully.
			LinpackSingular = -5, ///< linpack detected a computationally singular
								  ///< system using the symmetric indefinite
							  	  ///< factorization.
			ConvergeFail = -6, ///< the conjugate gradient iteration failed to
							   ///< converge in 30 iterations. the probable
							   ///< cause is an indefinite or near singular
							   ///< system. try using iflag=4. note that tol
							   ///< returns an estimate of the residual in
							   ///< the current conjugate gradient iteration.
	};

protected:
	vtkDbiharPatchFilter();
	~vtkDbiharPatchFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkDbiharPatchFilter(const vtkDbiharPatchFilter&); // Not implemented.
	void operator=(const vtkDbiharPatchFilter&); // Not implemented.

	double A, B, C, D;
	int MDim, NDim;
	int IFlag, OFlag;
	double Alpha, Beta;
	double Tol;
	int ITCG;

};

#endif
