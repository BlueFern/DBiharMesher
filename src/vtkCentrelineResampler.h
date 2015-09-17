#ifndef __vtkCentrelineResampler_h
#define __vtkCentrelineResampler_h

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

/**
 * This filter resamples lines in a centreline to ensure they all have an odd number of points
 * and are an even (and specified) distance apart.
 *
 * \param vtkPolyData A centreline.
 *
 * \param EdgeLength A number that specifies how far apart the resampled points should be.
 *
 * \return vtkPolyData with resampled points.
 */
class vtkCentrelineResampler : public vtkPolyDataAlgorithm
{

public:

	vtkTypeMacro(vtkCentrelineResampler, vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkCentrelineResampler *New();

	vtkSetMacro(EdgeLength, double);

	int GetNumbefOfBifurcations();

protected:

	vtkCentrelineResampler();
	~vtkCentrelineResampler() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

private:
	vtkCentrelineResampler(const vtkCentrelineResampler&); // Not implemented.
	void operator=(const vtkCentrelineResampler&); // Not implemented.

	double EdgeLength;

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
