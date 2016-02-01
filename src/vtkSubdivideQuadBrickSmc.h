#ifndef __vtkSubdivideQuadBrickSmc_h_
#define __vtkSubdivideQuadBrickSmc_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

/**
 * This filter subdivides a given quad into Columns x Rows smaller quads following a brick tessellation.
 * It uses vtkParametricSpline to interpolate points between the quads edges, converts
 * the newly created points into a structured grid and a geometry filter to build a mesh.
 *
 * The filter is specific to SMCs, gaps that are later filled by cells in neighbouring quads (following
 * the same subdivded pattern) are either on the top or the bottom of the input quad, depending on
 * the set value of the parameter Rotated. On the opposite side are 'short cells' (half sized), which are later stretched
 * into the gaps of the other circumferentially neighbouring quad.
 *
 * To later stretch the short cells across a second quad (relatively non-planar to the first) 6 points is required.
 * So that the cell can be modified in place, it too needs 6 points (vtk requirement). For consistency between this and
 * the similar filter for EC brick subdivision, two points of the short cell are duplicated and these points are on
 * the temporary edge (to be stretched later) This is important to note when connecting cells later in the pipeline.
 */
class vtkSubdivideQuadBrickSmc : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkSubdivideQuadBrickSmc,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkSubdivideQuadBrickSmc *New();

	vtkSetMacro(Columns, int);
	vtkSetMacro(Rows, int);
	vtkSetMacro(ShowProgress, bool);
	vtkSetMacro(Rotated, bool);

protected:
	vtkSubdivideQuadBrickSmc();
	~vtkSubdivideQuadBrickSmc() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	void copyPointsArray(double array1[], double array2[]);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkSubdivideQuadBrickSmc(const vtkSubdivideQuadBrickSmc&); // Not implemented.
	void operator=(const vtkSubdivideQuadBrickSmc&); // Not implemented.

	int Columns;
	int Rows;
	bool ShowProgress;
	bool Rotated;
};

#endif
