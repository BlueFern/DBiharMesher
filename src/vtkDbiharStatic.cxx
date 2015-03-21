#include <cmath>
#include <vtkMath.h>
#include <vtkIdList.h>
#include <vtkSmartPointer.h>

#include "vtkDbiharStatic.h"

const char *vtkDbiharStatic::RADII_VECTORS_ARR_NAME = {"radiiVectors"};
const char *vtkDbiharStatic::RADII_SCALARS_ARR_NAME = {"radiiScalars"};
const char *vtkDbiharStatic::DERIV_ARR_NAME = {"derivVectors"};
const char *vtkDbiharStatic::CELL_DATA_ARR_NAME = {"branchId"};

// The following constants are in micrometres.
const int vtkDbiharStatic::SMC_CIRC = 50;
const int vtkDbiharStatic::SMC_AXIAL = 5;
const int vtkDbiharStatic::EC_CIRC = 10;
const int vtkDbiharStatic::EC_AXIAL = 65;

void vtkDbiharStatic::DoubleCross1(const double v0[3], const double c0[3], const double v1[3], double c1[3])
{
	vtkMath::Cross(c0, v0, c1);
	vtkMath::Cross(v1, c1, c1);
}

/**
 * Find a global ID in a given IdList and return it's local position. Return -1 if not found.
 *
 * TODO: This method is redundant as there is equivalent functionality in vtkIdType vtkIdList::IsId(vtkIdType vtkid).
 *
 */
vtkIdType vtkDbiharStatic::GetPosition(vtkSmartPointer<vtkIdList> IdList, vtkIdType id)
{
	vtkIdType localId = 0;
	vtkIdType pointId = 0;
	vtkIdType maxLocalId = IdList->GetNumberOfIds() - 1;

	while(true)
	{
		if(localId > maxLocalId)
		{
			return -1;
		}
		pointId = IdList->GetId(localId);
		if(pointId == id)
		{
			break;
		}
		localId++;
	}
	return localId;
}
