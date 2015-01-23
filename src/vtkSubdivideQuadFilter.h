#ifndef __vtkSubdivideQuadFilter_h_
#define __vtkSubdivideQuadFilter_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

class vtkSubdivideQuadFilter : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkSubdivideQuadFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkSubdivideQuadFilter *New();

	vtkSetMacro(Columns, int);
	vtkSetMacro(Rows, int);


protected:
	vtkSubdivideQuadFilter();
	~vtkSubdivideQuadFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	void copyPointsArray(double array1[], double array2[]);

private:
	vtkSubdivideQuadFilter(const vtkSubdivideQuadFilter&); // Not implemented.
	void operator=(const vtkSubdivideQuadFilter&); // Not implemented.

	int Columns;
	int Rows;


};

#endif
