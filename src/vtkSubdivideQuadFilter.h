#ifndef __vtkSubdivideQuadFilter_h_
#define __vtkSubdivideQuadFilter_h_

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
class vtkSubdivideQuadFilter : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkSubdivideQuadFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkSubdivideQuadFilter *New();

	vtkSetMacro(Columns, int);
	vtkSetMacro(Rows, int);

protected:
	vtkSubdivideQuadFilter();
	~vtkSubdivideQuadFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	void copyPointsArray(double array1[], double array2[]);

private:
	vtkSubdivideQuadFilter(const vtkSubdivideQuadFilter&); // Not implemented.
	void operator=(const vtkSubdivideQuadFilter&); // Not implemented.

	int Columns;
	int Rows;


};

#endif
