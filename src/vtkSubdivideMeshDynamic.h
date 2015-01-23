#ifndef __vtkSubdivideMeshDynamic_h_
#define __vtkSubdivideMeshDynamic_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

class vtkSubdivideMeshDynamic : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkSubdivideMeshDynamic,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkSubdivideMeshDynamic *New();

	vtkSetMacro(Height, double);
	vtkSetMacro(Length, double);

protected:
	vtkSubdivideMeshDynamic();
	~vtkSubdivideMeshDynamic() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);


private:
	vtkSubdivideMeshDynamic(const vtkSubdivideMeshDynamic&); // Not implemented.
	void operator=(const vtkSubdivideMeshDynamic&); // Not implemented.

	double Height;
	double Length;


};

#endif
