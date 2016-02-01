#ifndef __vtkJoinSmcBrickMeshStatic_h_
#define __vtkJoinSmcBrickMeshStatic_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

/**
 * This filter accepts an SMC input mesh (that exists with gaps) and connects the cells.
 * It acts as a 'second pass' after each quad was subdivided with a brick tessellation
 * (which leaves holes at one end).
 *
 * The filter assumes the previous use of vtkSubdivideQuadBrickSmc, such that the top half (circumferentially) of
 * each branch was created with the parameter Rotated set to true, and the bottom half false. This difference
 * is necessary for a proper tessellation around the bifurcation (for a mesh that has one).
 *
 * It also leads to the circumferential edges having gaps - a problem particularly for non-flat
 * meshes. This is solved by stretching cells up and down from the lower and upper quad on this boundary,
 * respectively.
 *
 * Viewing the output from this filter directly will show some irregularities - cells that have been joined
 * so that they are non-plannar and have concave properties are not rendered correctly. OpenGL alone can't
 * handle non-convex polygons. This can be corrected by converting such polygons into two quads for example,
 * but should be done at a later stage and not in this filter as it would incorrectly increase the number of cells
 * The filter takes the number of axial and circumferential quads, the number of rows and columns each quad
 * was divided into, the number of branches, if the mesh is flat, and the mesh itself as parameters.
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
