
#ifndef __vtkSkipSegmentFilter_h
#define __vtkSkipSegmentFilter_h

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

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
