#ifndef __vtkJoinSmcBrickMeshStatic_h_
#define __vtkJoinSmcBrickMeshStatic_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>


class vtkJoinSmcBrickMesh : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkJoinSmcBrickMesh, vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkJoinSmcBrickMesh *New();

	vtkSetMacro(Columns, int);
	vtkSetMacro(Rows, int);
	vtkSetMacro(AxialQuads, int);
	vtkSetMacro(CircQuads, int);
	vtkSetMacro(Flat, bool);
	vtkSetMacro(Branches, int);

protected:
	vtkJoinSmcBrickMesh();
	~vtkJoinSmcBrickMesh() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkJoinSmcBrickMesh(const vtkJoinSmcBrickMesh&); // Not implemented.
	void operator=(const vtkJoinSmcBrickMesh&); // Not implemented.

	int Columns;
	int Rows;
	int AxialQuads;
	int CircQuads;
	bool Flat;
	int Branches;

};

#endif
