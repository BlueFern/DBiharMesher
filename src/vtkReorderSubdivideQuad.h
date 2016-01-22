#ifndef __vtkReorderSubdivideQuad_h_
#define __vtkReorderSubdivideQuad_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

/**
 * This filter subdivides a given quad into Columns x Rows smaller quads, where Columns and Rows are
 * parameters. It uses vtkParametricSpline to interpolate points between the quads edges, converts
 * the newly created points into a structured grid and a geometry filter to build a mesh.
 *
 * \param vtkPolydata A single quad - 4 points and one cell that is the quad itself.
 *
 * \param Columns Specifies the number of columns the original quad should be divided into.
 *
 * \param Rows Specifies the number of rows the original quad should be divided into.
 *
 * \return vtkPolyData containing Columns x Rows cells (as quads) and the new associated points.
 */
class vtkReorderSubdivideQuad : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkReorderSubdivideQuad,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkReorderSubdivideQuad *New();

	vtkSetMacro(Columns, int);
	vtkSetMacro(Rows, int);
	vtkSetMacro(Rotations, int);

protected:
	vtkReorderSubdivideQuad();
	~vtkReorderSubdivideQuad() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkReorderSubdivideQuad(const vtkReorderSubdivideQuad&); // Not implemented.
	void operator=(const vtkReorderSubdivideQuad&); // Not implemented.

	int Columns;
	int Rows;
	int Rotations;
};

#endif
