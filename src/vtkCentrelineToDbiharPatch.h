#ifndef __vtkCentrelineToDbiharPatch_h_
#define __vtkCentrelineToDbiharPatch_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPriorityQueue.h>

class vtkIdList;

class vtkCentrelineToDbiharPatch : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkCentrelineToDbiharPatch,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkCentrelineToDbiharPatch *New();

	vtkSetMacro(NumberOfRadialQuads, unsigned int);
	vtkSetMacro(SpineId, unsigned int);

protected:
	vtkCentrelineToDbiharPatch();
	~vtkCentrelineToDbiharPatch() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkCentrelineToDbiharPatch(const vtkCentrelineToDbiharPatch&); // Not implemented.
	void operator=(const vtkCentrelineToDbiharPatch&); // Not implemented.

	unsigned int NumberOfRadialQuads;
	unsigned int SpineId;
	double cEdgeScaling;
	double yEdgeScaling;

};

#endif
