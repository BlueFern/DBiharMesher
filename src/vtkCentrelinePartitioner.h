#ifndef __vtkCentrelinePartitioner_h_
#define __vtkCentrelinePartitioner_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

class vtkCentrelinePartitioner : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkCentrelinePartitioner,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkCentrelinePartitioner *New();
	static const int minEdgePoints;
	vtkSetMacro(PartitionLength, int);

protected:
	vtkCentrelinePartitioner();
	~vtkCentrelinePartitioner() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkCentrelinePartitioner(const vtkCentrelinePartitioner&); // Not implemented.
	void operator=(const vtkCentrelinePartitioner&); // Not implemented.

	void joinIdLists(vtkSmartPointer<vtkIdList> previous, vtkSmartPointer<vtkIdList> current,
					   vtkSmartPointer<vtkIdList> joined);
	void reverseIdList(vtkSmartPointer<vtkIdList> spine, vtkSmartPointer<vtkIdList> reversedSpine);

	int PartitionLength;
};

#endif
