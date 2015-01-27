#ifndef __vtkDbiharStatic_h_
#define __vtkDbiharStatic_h_

#include "vtkObject.h"

class vtkIdList;

class vtkDbiharStatic : public vtkObject {
public:

	static double AngleBetweenVectors(const double v1[3], const double v2[3]);
	static void DoubleCross1(const double v0[3], const double c0[3], const double v1[3], double c1[3]);
	static vtkIdType GetPosition(vtkSmartPointer<vtkIdList> IdList,vtkIdType id);

	static const char *RADII_VECTORS_ARR_NAME;
	static const char *RADII_SCALARS_ARR_NAME;
	static const char *DERIV_ARR_NAME;
	static const char *CELL_DATA_ARR_NAME;

	static const int SMC_CIRC;
	static const int SMC_AXIAL;
	static const int EC_CIRC;
	static const int EC_AXIAL;

};

#endif
