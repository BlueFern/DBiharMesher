#ifndef __vtkSubdivideQuadBrick_h_
#define __vtkSubdivideQuadBrick_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

class vtkSubdivideQuadBrick : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkSubdivideQuadBrick,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkSubdivideQuadBrick *New();

	vtkSetMacro(Columns, int);
	vtkSetMacro(Rows, int);
	vtkSetMacro(ShowProgress, bool);
	vtkSetMacro(Filled, bool);

protected:
	vtkSubdivideQuadBrick();
	~vtkSubdivideQuadBrick() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	void copyPointsArray(double array1[], double array2[]);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkSubdivideQuadBrick(const vtkSubdivideQuadBrick&); // Not implemented.
	void operator=(const vtkSubdivideQuadBrick&); // Not implemented.

	int Columns;
	int Rows;
	int CellType;
	bool ShowProgress;
	bool Filled;

};

#endif
