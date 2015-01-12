/*
 * Program: vtkPointsToMeshFilter.
 */

#ifndef __vtkPointsToMeshFilter_h
#define __vtkPointsToMeshFilter_h

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkUnsignedIntArray.h>

class vtkPointsToMeshFilter : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkPointsToMeshFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkPointsToMeshFilter *New();
	static const char *CELL_DATA_ARR_NAME;

	vtkSetObjectMacro(Dimensions, vtkUnsignedIntArray);

protected:
	vtkPointsToMeshFilter();
	~vtkPointsToMeshFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkPointsToMeshFilter(const vtkPointsToMeshFilter&); // Not implemented.
	void operator=(const vtkPointsToMeshFilter&); // Not implemented.

	vtkSmartPointer<vtkUnsignedIntArray> Dimensions;

};

#endif
