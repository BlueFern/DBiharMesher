#ifndef __vtkSkipSegmentFilter_h
#define __vtkSkipSegmentFilter_h

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

/**
 * This filter accepts a centreline and builds a skip segment mesh from a specified point
 * outwards. A skip segment is a section of the mesh that is generated purely for computational
 * fluid dynamics; it follows an imaginary/newly created section of the centreline. The direction
 * and spacing between points for this generated centreline is determined by point Id and the point
 * after or before it (for an Inlet or Outlet respectively).
 *
 * \param vtkPolyData Centreline data.
 *
 * \param SkipSize Determines how many rings of points the skip segment will contain.
 *
 * \param Inlet Boolean flag that specifies if the given point Id is positioned at an
 * inlet or not.
 *
 * \param Outlet Boolean flag that specifies if the given point Id is positioned at an
 * outlet or not.
 *
 * \param PointId Id of the point the start of the skip segment should be centred around.
 *
 * \returns vtkPolyData that contains points around the skip segment, a triangular mesh
 * around the surface of this segment and one line that that is the final ring of points
 * (which is useful later for building caps).
 */
class vtkSkipSegmentFilter : public vtkPolyDataAlgorithm
{
public:
	vtkTypeMacro(vtkSkipSegmentFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkSkipSegmentFilter *New();

	vtkSetMacro(SkipSize, int);
	vtkSetMacro(NumberOfRadialQuads, int);
	vtkSetMacro(Inlet, bool);
	vtkSetMacro(Outlet, bool);
	vtkSetMacro(PointId, vtkIdType);

protected:
	vtkSkipSegmentFilter();
	~vtkSkipSegmentFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkSkipSegmentFilter(const vtkSkipSegmentFilter&); // Not implemented.
	void operator=(const vtkSkipSegmentFilter&); // Not implemented.
	int NumberOfRadialQuads;
	int SkipSize;
	bool Inlet;
	bool Outlet;
	vtkIdType PointId;

};

#endif
