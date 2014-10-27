/*
 * Program: vtkDbiharPatchFilter.
 */

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

/*
 * vtkDbiharPatchFilter is a wrapper for NetLib's biharmonic PDF equation solver.
 * See <a href="http://www.netlib.org/bihar/dbihar.f">Dbihar on NetLib</a>.
 * The wrapper is written as a VTK-style filter.
 */
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
	vtkSetMacro(MQuads,int);
	vtkSetMacro(NQuads,int);
	vtkSetMacro(Tol,double);
	vtkSetMacro(IFlag,int);
	vtkSetMacro(Alpha,double);
	vtkSetMacro(Beta,double);

	static const char *DERIV_ARR_NAME;

protected:
	vtkDbiharPatchFilter();
	~vtkDbiharPatchFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkDbiharPatchFilter(const vtkDbiharPatchFilter&); // Not implemented.
	void operator=(const vtkDbiharPatchFilter&); // Not implemented.

	double A, B, C, D;
	int MDim, NDim;
	int MQuads, NQuads;
	int IFlag, OFlag;
	double Alpha, Beta;
	double Tol;
	int ITCG;

	void CheckError();
};

#endif
