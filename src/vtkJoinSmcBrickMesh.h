#ifndef __vtkJoinSmcBrickMeshStatic_h_
#define __vtkJoinSmcBrickMeshStatic_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

/**
 * This filter accepts an SMC input mesh (that exists with gaps) and connects the cells.
 * It acts as a 'second pass' after each quad was subdivided with a brick tessellation
 * (which leaves holes at one end).
 *
 * The filter assumes the previous use of vtkReorderSubdivideQuad, that the top half of each
 * branch was created with a rotation value of 3 (90 degrees counter-clockwise), and the bottom
 * half with 1 (90 degrees clockwise). This difference is necessary for a proper tessellation around
 * the bifurcation (for a mesh that has one).
 *
 * This difference leads to the circumferential edges having gaps - a problem particularly for a non-flat
 * mesh. This is solved by stretching cells up and down from the lower and upper quad on this boundary,
 * respectively.
 *
 * Viewing this output shows some irregularities - cells that have been joined so that they are non-plannar
 * and have concave properties are not rendered correctly. OpenGL alone can't handle non-convex polygons.
 * For future considerations, this is not addressed in this filter. At a later point this mesh will be put through
 * a triangle filter which solves this issue.
 *
 * The filter takes the number of axial and circumferential quads, the number of rows and columns each quad
 * was divided into, the number of branches, if the mesh is flat, and the mesh itself as parameters.
 *
 * The output is a fully joined mesh.
 *
 */
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
