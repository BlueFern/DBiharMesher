#ifndef __vtkSubdivideQuadBrick_h_
#define __vtkSubdivideQuadBrick_h_

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
class vtkSubdivideQuadBrick : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkSubdivideQuadBrick,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkSubdivideQuadBrick *New();

	vtkSetMacro(Columns, int);
	vtkSetMacro(Rows, int);
	vtkSetMacro(ShowProgress, bool);
	vtkSetMacro(CellType, int);

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

};

#endif
