#ifndef __vtkJoinEcBrickMeshStatic_h_
#define __vtkJoinEcBrickMeshStatic_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

/**
 * This filter accepts an EC input mesh (that exists with gaps) and connects the cells.
 * It acts as a 'second pass' after each quad was subdivided with a brick tessellation
 * (which leaves holes at one end).
 *
 * For a correct tessellation, the first of three branches (for such a mesh) uses a rotation value of 0
 * (no rotation), while the other two use 2 (180 degrees clockwise) in the vtkReorderSubdivdeQuad filter.
 * The input mesh for this filter is assumed to have followed this restriction.
 *
 * Viewing the output from this filter directly will show some irregularities - cells that have been joined
 * so that they are non-plannar and have concave properties are not rendered correctly. OpenGL alone can't
 * handle non-convex polygons. For future considerations, this is not addressed in this filter. At a later point
 * this mesh will be put through a triangle filter which solves this issue.
 *
 * The filter takes the number of axial and circumferential quads, the number of rows and columns each quad
 * was divided into, the number of branches, if the mesh is flat, and the mesh itself as parameters.
 *
 * The output is a fully joined mesh.
 *
 */
class vtkJoinEcBrickMesh : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkJoinEcBrickMesh, vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkJoinEcBrickMesh *New();

	vtkSetMacro(Columns, int);
	vtkSetMacro(Rows, int);
	vtkSetMacro(AxialQuads, int);
	vtkSetMacro(CircQuads, int);
	vtkSetMacro(Branches, int);

protected:
	vtkJoinEcBrickMesh();
	~vtkJoinEcBrickMesh() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkJoinEcBrickMesh(const vtkJoinEcBrickMesh&); // Not implemented.
	void operator=(const vtkJoinEcBrickMesh&); // Not implemented.

	int Columns;
	int Rows;
	int AxialQuads;
	int CircQuads;
	int Branches;

};

#endif
