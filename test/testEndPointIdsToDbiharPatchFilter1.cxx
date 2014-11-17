
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

#include <vtkXMLStructuredGridReader.h>
#include <vtkStructuredGridAppend.h>

#include "showPolyData.h"

#include "vtkCentrelineData.h"
#include "vtkScalarRadiiToVectorsFilter.h"
#include "vtkEndPointIdsToDbiharPatchFilter.h"

#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/SyntheticBifurcation_1.vtk").c_str());
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

#if 1
	endPointIdsList->InsertNextId(40);
	endPointIdsList->InsertNextId(100);
	endPointIdsList->InsertNextId(176);
#else
	endPointIdsList->InsertNextId(0);
	endPointIdsList->InsertNextId(70);
#endif

	const double unitsConversionFactor = 1.0e-3;
	vtkSmartPointer<vtkDoubleArray> radiiArray = vtkDoubleArray::SafeDownCast(resampledVesselCentreline->GetPointData()->GetArray(vtkCentrelineData::RADII_ARR_NAME));
	double R = radiiArray->GetValue(endPointIdsList->GetId(0));
	R *= unitsConversionFactor; // R in m.

	double C = 2 * vtkMath::Pi() * R;

	const double ECLength = 65e-6; // m.
	const double SMCLength = 50e-6; // m.
	const unsigned int ECMultiple = 4;
	const unsigned int SMCMultiple = 4;

	double tmpDoubleVal = (C / 2.0) / (SMCLength * SMCMultiple);
	int tmpIntVal = vtkMath::Round(tmpDoubleVal);

	int numberOfRadialQuads = tmpIntVal;
	// Must be even.
	if((numberOfRadialQuads & 1) == 1)
	{
		if(tmpDoubleVal - (double)tmpIntVal > 0)
		{
			// Increment.
			numberOfRadialQuads++;
		}
		else
		{
			// Decrement.
			numberOfRadialQuads--;
		}
	}

	vtkSmartPointer<vtkEndPointIdsToDbiharPatchFilter> idListToDbiharPatchFilter = vtkSmartPointer<vtkEndPointIdsToDbiharPatchFilter>::New();
	idListToDbiharPatchFilter->SetInputData(resampledVesselCentrelineWithRadii);
	idListToDbiharPatchFilter->SetNumberOfRadialQuads(numberOfRadialQuads);
	idListToDbiharPatchFilter->SetEndPointIdsList(endPointIdsList);
	idListToDbiharPatchFilter->Update();

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}

