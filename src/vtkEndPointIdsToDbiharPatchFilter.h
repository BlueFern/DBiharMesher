/**
 * Program: vtkEndPointIdsToDbiharPatchFilter.
 */

#ifndef __vtkEndPointIdsToDbiharPatchFilter_h
#define __vtkEndPointIdsToDbiharPatchFilter_h

#include <map>

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

class vtkEndPointIdsToDbiharPatchFilter: public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkEndPointIdsToDbiharPatchFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkEndPointIdsToDbiharPatchFilter *New();

	vtkSetObjectMacro(EndPointIdsList, vtkIdList);

protected:
	vtkEndPointIdsToDbiharPatchFilter();
	~vtkEndPointIdsToDbiharPatchFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkEndPointIdsToDbiharPatchFilter(const vtkEndPointIdsToDbiharPatchFilter&); // Not implemented.
	void operator=(const vtkEndPointIdsToDbiharPatchFilter&); // Not implemented.

	vtkPolyData* input;
	std::map<vtkIdType, std::vector<vtkIdType> > treeInfo;
	vtkSmartPointer<vtkIdList> EndPointIdsList;

	std::vector<std::vector<vtkIdType> > spineIds;
};

#endif
