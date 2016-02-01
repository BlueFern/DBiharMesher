#ifndef __vtkSubdivideMeshBrickStatic_h_
#define __vtkSubdivideMeshBrickStatic_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

/**
 * This filter subdivides all quads in the input into ECs or SMCs depending on
 * the CellType parameter. It then iterates over the generated mesh (that will exist
 * with gaps due to the brick tessellation) and connect cells between quads.
 *
 * The output is then either a full EC or SMC mesh.
 */
class vtkSubdivideMeshBrick : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkSubdivideMeshBrick, vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkSubdivideMeshBrick *New();

	vtkSetMacro(Columns, int);
	vtkSetMacro(Rows, int);
	vtkSetMacro(AxialQuads, int);
	vtkSetMacro(CircQuads, int);
	vtkSetMacro(Flat, bool);
	vtkSetMacro(CellType, int);
	vtkSetMacro(Branches, int);

protected:
	vtkSubdivideMeshBrick();
	~vtkSubdivideMeshBrick() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkSubdivideMeshBrick(const vtkSubdivideMeshBrick&); // Not implemented.
	void operator=(const vtkSubdivideMeshBrick&); // Not implemented.

	int Columns;
	int Rows;
	int AxialQuads;
	int CircQuads;
	int CellType;
	int Branches;

	bool Flat;

};

#endif
