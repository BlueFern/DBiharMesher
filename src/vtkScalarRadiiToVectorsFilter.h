/**
 * Program: vtkScalarRadiiToVectorsFilter.
 */

#ifndef __vtkScalarRadiiToVectorsFilter_h
#define __vtkScalarRadiiToVectorsFilter_h

#include "vtkPolyDataAlgorithm.h"

class vtkScalarRadiiToVectorsFilter: public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkScalarRadiiToVectorsFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkScalarRadiiToVectorsFilter *New();

protected:
	vtkScalarRadiiToVectorsFilter();
	~vtkScalarRadiiToVectorsFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkScalarRadiiToVectorsFilter(const vtkScalarRadiiToVectorsFilter&); // Not implemented.
	void operator=(const vtkScalarRadiiToVectorsFilter&); // Not implemented.

};

#endif
