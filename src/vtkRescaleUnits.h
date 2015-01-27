#ifndef __vtkRescaleUnits_h_
#define __vtkRescaleUnits_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

class vtkRescaleUnits : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkRescaleUnits,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkRescaleUnits *New();

	vtkSetMacro(Scale, int);

protected:
	vtkRescaleUnits();
	~vtkRescaleUnits() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkRescaleUnits(const vtkRescaleUnits&); // Not implemented.
	void operator=(const vtkRescaleUnits&); // Not implemented.

	int Scale;
};

#endif
