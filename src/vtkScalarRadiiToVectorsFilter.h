#ifndef __vtkScalarRadiiToVectorsFilter_h
#define __vtkScalarRadiiToVectorsFilter_h

#include <map>

#include <vtkVector.h>
#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

/**
 * From the radiiScalars point data this filter creates radiiVectors that are vectors from each point
 * in the centreline out to the surface of the vessel.
 *
 * \param vtkPolyData A centreline.
 *
 * \param angleTolerance An optional parameter that specifies the tolerance between angles.
 *
 * \return vtkPolyData with an extra array of point data, radiiVectors.
 */
class vtkScalarRadiiToVectorsFilter: public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkScalarRadiiToVectorsFilter,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkScalarRadiiToVectorsFilter *New();

protected:
	vtkScalarRadiiToVectorsFilter();
	~vtkScalarRadiiToVectorsFilter() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkScalarRadiiToVectorsFilter(const vtkScalarRadiiToVectorsFilter&); // Not implemented.
	void operator=(const vtkScalarRadiiToVectorsFilter&); // Not implemented.

	double angleTolerance;

	vtkPolyData* inputPointerCopy;

	std::map<vtkIdType, std::vector<vtkIdType> > treeInfo;
	std::map<vtkIdType, vtkVector3d> avrgVectors;

	void GetDirectionVector(vtkIdType lineId, vtkIdType pointId, double *vector);

	vtkSmartPointer<vtkIdList> GetLineIds(vtkIdType lineId);
};

#endif
