/**
 * Program: vtkScalarRadiiToVectorsFilter.
 */

#ifndef __vtkScalarRadiiToVectorsFilter_h
#define __vtkScalarRadiiToVectorsFilter_h

#include <map>

#include <vtkVector.h>
#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

// TODO: Declare this within the class.
typedef enum {d_start = 0, d_end = 1} LocationType;

class vtkScalarRadiiToVectorsFilter: public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkScalarRadiiToVectorsFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkScalarRadiiToVectorsFilter *New();

	static const char *RADII_ARR_NAME;

protected:
	vtkScalarRadiiToVectorsFilter();
	~vtkScalarRadiiToVectorsFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkScalarRadiiToVectorsFilter(const vtkScalarRadiiToVectorsFilter&); // Not implemented.
	void operator=(const vtkScalarRadiiToVectorsFilter&); // Not implemented.

	double angleTolerance;

	vtkPolyData* input;

	std::map<vtkIdType, std::vector<vtkIdType> > treeInfo;
	std::map<vtkIdType, vtkVector3d> avrgVectors;

	//void GetDirectionVector(vtkIdType lineId, LocationType location, double *vector);
	void GetDirectionVector(vtkIdType lineId, vtkIdType pointId, double *vector);

	//vtkIdType GetGlobalPointId(vtkIdType lineId, LocationType location);
	//vtkIdType GetGlobalPointId(vtkIdType lineId, vtkIdType pointId);

	vtkSmartPointer<vtkIdList> GetLineIds(vtkIdType lineId);

};

#endif
