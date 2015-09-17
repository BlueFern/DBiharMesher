#include <map>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkUnsignedIntArray.h>
#include <vtkCellArray.h>

#include <vtkPointData.h>
#include <vtkXMLStructuredGridReader.h>

#include "vtkDbiharStatic.h"
#include "vtkCentrelineToDbiharPatch.h"
#include "vtkCentrelinePartitioner.h"
#include "vtkScalarRadiiToVectorsFilter.h"

#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/testCentrelineToDbiharPatch0.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(vesselCentreline);
	scalarRadiiToVectorsFilter->Update();

	vtkSmartPointer<vtkCentrelinePartitioner> centrelinePartitioner = vtkSmartPointer<vtkCentrelinePartitioner>::New();
	centrelinePartitioner->SetInputData(scalarRadiiToVectorsFilter->GetOutput());
	centrelinePartitioner->SetPartitionLength(50);
	vtkSmartPointer<vtkIdList> EndPoints = vtkSmartPointer<vtkIdList>::New();
	EndPoints->InsertNextId(20); // Remove the curve in inlet, causes folding in mesh otherwise.

	centrelinePartitioner->SetEndPoints(EndPoints);
	centrelinePartitioner->Update();

	vtkSmartPointer<vtkCentrelineToDbiharPatch> dbiharPatchFilter = vtkSmartPointer<vtkCentrelineToDbiharPatch>::New();

	dbiharPatchFilter->SetShowProgress(true);
	dbiharPatchFilter->SetInputData(centrelinePartitioner->GetOutput());
	dbiharPatchFilter->SetNumberOfRadialQuads(28);
	dbiharPatchFilter->SetSpineId(0);
	dbiharPatchFilter->SetArchDerivScale(3.2);
	dbiharPatchFilter->SetEdgeDerivScale(4.0);
	dbiharPatchFilter->Update();

	vtkDbiharStatic::WritePolyData(dbiharPatchFilter->GetOutput(), "centrelineToDbiharTest.vtp");

	return EXIT_SUCCESS;
}

