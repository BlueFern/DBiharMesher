/*
 * Program: vtkCentrelineData.
 */

#ifndef __vtkCentrelineData_h
#define __vtkCentrelineData_h

#include <vtkPolyData.h>

/*
 * vtkCentrelineData is a wrapper for vtkPolyData for extracting segments of data from centrelines.
 */
class vtkCentrelineData : public vtkPolyData
{
public:
	vtkTypeMacro(vtkCentrelineData,vtkPolyData);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkCentrelineData *New();

	int GetNumVessels();
	int GetNumBifurcations();

protected:
	vtkCentrelineData();
	~vtkCentrelineData() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkCentrelineData(const vtkCentrelineData&); // Not implemented.
	void operator=(const vtkCentrelineData&); // Not implemented.

};

#endif
