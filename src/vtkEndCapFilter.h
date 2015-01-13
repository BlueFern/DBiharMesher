#ifndef __vtkEndCapFilter_h
#define __vtkEndCapFilter_h

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

class vtkEndCapFilter : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkEndCapFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkEndCapFilter *New();

protected:
	vtkEndCapFilter();
	~vtkEndCapFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkEndCapFilter(const vtkEndCapFilter&); // Not implemented.
	void operator=(const vtkEndCapFilter&); // Not implemented.
};

#endif
