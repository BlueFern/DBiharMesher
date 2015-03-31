#ifndef __vtkSubdivideMeshStatic_h_
#define __vtkSubdivideMeshStatic_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

/**
 * This filter subdivides all quads in the input to cells roughly of Height x Length size.
 * For every quad in the input vtkSubdivideQuadFilter is used, where the number of columns
 * and rows are calculated depending on how many smaller cells can fit into each original
 * quad. This would mean, for example, at narrower sections of the vessel there will be fewer
 * new quads created because the original quads are smaller.
 *
 * \param vtkPolyData A quad mesh. This filter will subdivide all existing cells
 * in this input.
 *
 * \param Height The circumferential length of the new quads.
 *
 * \param Length The axial length of the new quads.
 *
 * \return vtkPolyData with cells of roughly Height x Length size.
 */
class vtkSubdivideMeshStatic : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkSubdivideMeshStatic, vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkSubdivideMeshStatic *New();

	vtkSetMacro(Columns, double);
	vtkSetMacro(Rows, double);

protected:
	vtkSubdivideMeshStatic();
	~vtkSubdivideMeshStatic() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkSubdivideMeshStatic(const vtkSubdivideMeshStatic&); // Not implemented.
	void operator=(const vtkSubdivideMeshStatic&); // Not implemented.

	double Columns;
	double Rows;


};

#endif
