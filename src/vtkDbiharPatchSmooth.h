#ifndef __vtkDbiharPatchSmooth_h_
#define __vtkDbiharPatchSmooth_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

class vtkIdList;


class vtkDbiharPatchSmooth : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkDbiharPatchSmooth,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkDbiharPatchSmooth *New();

	vtkSetMacro(NumRadialQuads, int);

protected:
	vtkDbiharPatchSmooth();
	~vtkDbiharPatchSmooth() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);


private:
	vtkDbiharPatchSmooth(const vtkDbiharPatchSmooth&); // Not implemented.
	void operator=(const vtkDbiharPatchSmooth&); // Not implemented.

	int NumRadialQuads;


};

#endif
