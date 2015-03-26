#ifndef __vtkCentrelineToDbiharPatch_h_
#define __vtkCentrelineToDbiharPatch_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPriorityQueue.h>

class vtkIdList;

/**
 * This filter collects the external points around a specified line section (the centreline
 * has assumed to have been previously partitioned using vtkCentrelinePartitioner). This collection
 * is sent to the Dbihar filter to create the internal points.
 *
 * Both Dirichlet and Neumann boundary conditions are stored in a vtkPolyData object.
 * The Dirichlet boundary values are stored as the ordered set of points.
 * The Neumann boundary values are stored as vectors associated with each corresponding point.
 * For a non-bifurcating surface mesh there are two list of point ids, one per input patch for the Dbihar filter.
 * For a bifurcating surface mesh there are three lists of point ids, one per input patch for the Dbihar filter.
 * The code walks the given point id list starting at the first point. For this centreline point id an arch consisting of points obtained through radius rotation around this centreline point is inserted into the Dirichlet boundary condition point set.
 * Then the centreline point id list is traversed forward, where for each point id and the corresponding radius vector a point is added into the Dirichlet boundary condition point set.
 * When the last centreline point id is reached, an arch consisting of points obtained through radius rotation around this centreline point is inserted into the Dirichlet boundary condition point set.
 * Then the centreline point id list is traversed backward, where for each point id and the corresponding radius vector a point is added into the Dirichlet boundary condition point set.
 * The Neumann boundary conditions are inserted for every point during the centreline point id list traversal. For the arches the direction of the derivatives is parallel to the direction of the centreline at the start and end points of the centreline segment.
 * For the other two edges of the patch the derivatives are oriented perpendicular to both the vector radius and the corresponding centreline edge. The magnitudes of the derivative vectors are proportional to the radius of the centreline at the given point.
 *
 * \param vtkPolyData Already partitioned centreline data.
 *
 * \param CellId A particular cell to build the boundary around.
 *
 * \return vtkPolyData. The points are the ones returned from vtkDbiharPatchFilter which is
 *  internally called after the boundary is found.
 */

class vtkCentrelineToDbiharPatch : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkCentrelineToDbiharPatch,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkCentrelineToDbiharPatch *New();

	vtkSetMacro(ArchDerivScale, double);
	vtkSetMacro(EdgeDerivScale, double);
	vtkSetMacro(NumberOfRadialQuads, unsigned int);
	vtkSetMacro(SpineId, unsigned int);
	vtkSetMacro(BifurcationId, unsigned int);
	vtkSetMacro(ShowProgress, bool);

protected:
	vtkCentrelineToDbiharPatch();
	~vtkCentrelineToDbiharPatch() {};

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkCentrelineToDbiharPatch(const vtkCentrelineToDbiharPatch&); // Not implemented.
	void operator=(const vtkCentrelineToDbiharPatch&); // Not implemented.

	unsigned int NumberOfRadialQuads;
	unsigned int SpineId;
	double ArchDerivScale;
	double EdgeDerivScale;
	unsigned int BifurcationId;
	bool ShowProgress;

};

#endif
