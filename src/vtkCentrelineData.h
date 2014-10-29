/*
 * Program: vtkCentrelineData.
 */

#ifndef __vtkCentrelineData_h
#define __vtkCentrelineData_h

class vtkPolyData;

// TODO: This class probably is better named as vtkCentrelineResampler.
// TODO: This class probably is better implemented as a filter.

/*
 * vtkCentrelineData is a wrapper for vtkPolyData for extracting segments of data from centrelines.
 */
class vtkCentrelineData : public vtkObject
{
public:
	vtkTypeMacro(vtkCentrelineData,vtkObject);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkCentrelineData *New();

	void SetCentrelineData(vtkPolyData *centrelineData);
	vtkPolyData *GetOutput();

	int GetNumbefOfBifurcations();

protected:
	vtkCentrelineData();
	~vtkCentrelineData() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkCentrelineData(const vtkCentrelineData&); // Not implemented.
	void operator=(const vtkCentrelineData&); // Not implemented.

	vtkSmartPointer<vtkPolyData> polyData;

	const double unitsConversionFactor = 1.0e-3;
	const double ECLength = 65e-6;
	const double SMCLength = 50e-6;
	const unsigned int ECMultiple = 4;
};

#endif
