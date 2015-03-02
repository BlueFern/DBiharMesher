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

protected:
	vtkCentrelineToDbiharPatch();
	~vtkCentrelineToDbiharPatch() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkCentrelineToDbiharPatch(const vtkCentrelineToDbiharPatch&); // Not implemented.
	void operator=(const vtkCentrelineToDbiharPatch&); // Not implemented.

	unsigned int NumberOfRadialQuads;
	unsigned int SpineId;
	double ArchDerivScale;
	double EdgeDerivScale;

};

#endif
