#ifndef __vtkPointsToMeshFilter_h
#define __vtkPointsToMeshFilter_h

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkUnsignedIntArray.h>

/**
 * This filter builds a quadrilateral mesh using points specified in input vtkPolyData.
 *
 * A vtkIntArray is created and used as cell data for the quads. It serves to
 * reference what trunk a particular cell belongs to; so for each quad created a value
 * is put into the array (0 if the cell belongs to the first trunk, 1 for the next, etc.).
 *
 * \param vtkPolyData This contains no cells or attribute data, only points that form
 * a 3D model of a straight segment or bifurcation. The order of the points is important;
 * they follow the spines around a bifurcation or up then down a straight segment in half loops/circles.
 *
 * \param Dimensions An array that gives details about the particular input vtkPolyData. The first
 * element specifies the number of interior quads in each half ring (so is one less than
 * the number of points in it). The next one or more elements are the number of rows of
 * quads that will fit in a branch (so one less than the number of half rings in it).
 * (There would therefore be Dimensions->GetValue(0) * Dimensions->GetValue(b) quads on
 * branch 'b'.)
 *
 * \return vtkPolyData is returned with the quads as cells and most of the points from the input.
 * Duplicate points between connecting halves are omitted, but the points between connected
 * branches are included for each branch it is associated with. This association is actually
 * fictitious, as there are no cells or point data, but is nonetheless made for simplicity.
 * The order of these points is also different. Instead of half rings the points move entirely
 * around the model. If we're not dealing with straight segments, points start at the bifurcation
 * and move out to the end of the branch.
 */
class vtkPointsToMeshFilter : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkPointsToMeshFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkPointsToMeshFilter *New();

	vtkSetObjectMacro(Dimensions, vtkUnsignedIntArray);
	vtkSetMacro(ShowProgress, bool);

protected:
	vtkPointsToMeshFilter();
	~vtkPointsToMeshFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkPointsToMeshFilter(const vtkPointsToMeshFilter&); // Not implemented.
	void operator=(const vtkPointsToMeshFilter&); // Not implemented.

	vtkSmartPointer<vtkUnsignedIntArray> Dimensions;
	bool ShowProgress;

};

#endif
