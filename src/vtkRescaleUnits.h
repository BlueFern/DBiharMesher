#ifndef __vtkRescaleUnits_h_
#define __vtkRescaleUnits_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

/**
 * Re-scales a data set by a specified amount. This includes both the points
 * and point data.
 *
 * \param vtkPolyData Centreline data.
 *
 * \param Scale Apecifies the amount to scale the input.
 *
 * \return vtkPolyData exactly the same as the input, only scaled by a given amount.
 *
 */
class vtkRescaleUnits : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkRescaleUnits,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkRescaleUnits *New();

	vtkSetMacro(Scale, int);

protected:
	vtkRescaleUnits();
	~vtkRescaleUnits() {};

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkRescaleUnits(const vtkRescaleUnits&); // Not implemented.
	void operator=(const vtkRescaleUnits&); // Not implemented.

	int Scale;
};

#endif
