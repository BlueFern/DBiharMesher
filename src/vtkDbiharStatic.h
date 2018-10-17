#ifndef __vtkDbiharStatic_h_
#define __vtkDbiharStatic_h_

#include "vtkObject.h"

class vtkIdList;
class vtkPolyData;
class vtkStructuredGrid;

/**
 * This is a class for storing string and numeric constants as well as shared static methods.
 */
class vtkDbiharStatic : public vtkObject {
public:

	static const char *RADII_VECTORS_ARR_NAME;
	static const char *RADII_SCALARS_ARR_NAME;
	static const char *DERIV_ARR_NAME;
	static const char *BRANCH_ID_ARR_NAME;
	static const char *GRID_COORDS_ARR_NAME;

	static const double SMC_CIRC;
	static const double SMC_AXIAL;
	static const double EC_CIRC;
	static const double EC_AXIAL;

	/**
	 * Get formatted date/time stamp as string.
	 */
	static std::string GetTimeStamp();

	/**
	 * Write polydata in XML format.
	 */
	static void WritePolyData(vtkPolyData *input, std::string fileName);

	/**
	 * Write polydata in STL format.
	 */
	static void WriteStlData(vtkPolyData *input, std::string fileName);

	/**
	 * Show input and output polydata.
	 */
	static void ShowPolyDataWithGrid(vtkPolyData *input, vtkStructuredGrid *output, double derivateScaling = 0.1);

	/**
	 * Show polydata with vectors.
	 */
	static void ShowPolyData(vtkPolyData *input, double vectorScaling = 1.0);

	static void DoubleCross(const double v0[3], const double c0[3], const double v1[3], double c1[3]);

	static void PrintDataArray(vtkDataArray *dataArray);

};

#endif
