
#include <map>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectReader.h>

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

#include "vtkCentrelineData.h"
#include "vtkScalarRadiiToVectorsFilter.h"
#include "showPolyData.h"

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

	//std::cout << "Number of input points: " << vesselCentreline->GetNumberOfPoints() << std::endl;
	//std::cout << "Number of input lines: " << vesselCentreline->GetNumberOfLines() << std::endl;

	//std::cout << "Number of output points: " << resampledVesselCentreline->GetNumberOfPoints() << std::endl;
	//std::cout << "Number of output lines: " << resampledVesselCentreline->GetNumberOfLines() << std::endl;

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(resampledVesselCentreline);
	scalarRadiiToVectorsFilter->Update();

	vtkPolyData *resampledVesselCentrelineWithRadii = scalarRadiiToVectorsFilter->GetOutput();
	showPolyData1(resampledVesselCentrelineWithRadii, 1.0);

#if 1
	vtkSmartPointer<vtkXMLPolyDataWriter> tmpWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	tmpWriter->SetInputData(resampledVesselCentrelineWithRadii);
	tmpWriter->SetFileName("resampledCentrelineWithRadii.vtp");
	tmpWriter->Write();

	vtkSmartPointer<vtkGenericDataObjectWriter> writer = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName("resampledCentrelineWithRadii.vtk");
	writer->SetInputData(resampledVesselCentrelineWithRadii);
	writer->Write();
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}

