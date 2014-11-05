
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

#include <vtkCellArray.h>
#include <vtkPointData.h>

#include <vtkXMLPolyDataWriter.h>
#include <vtkGenericDataObjectWriter.h>

#include "showPolyData.h"

#include "vtkCentrelineData.h"
#include "vtkScalarRadiiToVectorsFilter.h"
#include "vtkEndPointIdsToDbiharPatchFilter.h"

#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/227A_Centreline.vtk").c_str());
	//vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/721A_Centreline.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	vtkSmartPointer<vtkCentrelineData> centrelineSegmentSource = vtkSmartPointer<vtkCentrelineData>::New();
	centrelineSegmentSource->SetCentrelineData(vesselCentreline);

	vtkPolyData *resampledVesselCentreline = centrelineSegmentSource->GetOutput();

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(resampledVesselCentreline);
	scalarRadiiToVectorsFilter->Update();

	vtkPolyData *resampledVesselCentrelineWithRadii = scalarRadiiToVectorsFilter->GetOutput();

	vtkSmartPointer<vtkIdList> endPointIdsList = vtkSmartPointer<vtkIdList>::New();
	endPointIdsList->InsertNextId(21);
	endPointIdsList->InsertNextId(79);
	endPointIdsList->InsertNextId(948);

	const double unitsConversionFactor = 1.0e-3;
	vtkSmartPointer<vtkDoubleArray> radiiArray = vtkDoubleArray::SafeDownCast(resampledVesselCentreline->GetPointData()->GetArray(vtkCentrelineData::RADII_ARR_NAME));
	double R = radiiArray->GetValue(39);
	R *= unitsConversionFactor; // R in m.

	double C = 2 * vtkMath::Pi() * R;

	const double ECLength = 65e-6; // m.
	const double SMCLength = 50e-6; // m.
	const unsigned int ECMultiple = 4;
	const unsigned int SMCMultiple = 4;

	unsigned int numberOfRadialQuads = (C / 2.0) / (SMCLength * SMCMultiple);

	vtkSmartPointer<vtkEndPointIdsToDbiharPatchFilter> idListToDbiharPatchFilter = vtkSmartPointer<vtkEndPointIdsToDbiharPatchFilter>::New();
	idListToDbiharPatchFilter->SetInputData(resampledVesselCentrelineWithRadii);
	idListToDbiharPatchFilter->SetNumberOfRadialQuads(numberOfRadialQuads);
	idListToDbiharPatchFilter->SetEndPointIdsList(endPointIdsList);
	idListToDbiharPatchFilter->Update();

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}

