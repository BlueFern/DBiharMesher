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

	static const char *RADII_ARR_NAME;

	void SetCentrelineData(vtkPolyData *centrelineData);
	vtkSetMacro(EdgeLength, double);
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

	double EdgeLength = 0;

#if 0
	// TODO: unitsConversionFactor should be defined through a setter method, if at all. Perhaps it is best this class should be agnostic of these constants.
	const double unitsConversionFactor = 1.0e-3;
	// TODO: EC_Multiple and SMC_Multiple should be defined through a setter method, if at all. Perhaps it is best this class should be agnostic of these constants.
	const unsigned int ECMultiple = 4;
	const unsigned int SMCMultiple = 4;
	// TODO: These constants should be available globally. Perhaps it is best this class should be agnostic of these constants.
	const double ECLength = 65e-6;
	const double SMCLength = 50e-6;
#endif

};

#endif
