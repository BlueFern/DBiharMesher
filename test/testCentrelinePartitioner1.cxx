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
#include <vtkPolyDataWriter.h>

#include <vtkCellArray.h>
#include <vtkCentrelinePartitioner.h>
#include <vtkPointData.h>

#include <vtkXMLPolyDataWriter.h>
#include <vtkGenericDataObjectWriter.h>

#include <vtkXMLStructuredGridReader.h>

#include "showPolyData.h"

#include "vtkCentrelineData.h"
#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/227A_CentrelineResampled_4ECs.vtk").c_str());
	//vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/721A_CentrelineResampled_4ECs.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	//vtkSmartPointer<vtkCentrelineData> centrelineSegmentSource = vtkSmartPointer<vtkCentrelineData>::New();
	//centrelineSegmentSource->SetCentrelineData(vesselCentreline);

	//vtkPolyData *resampledVesselCentreline = centrelineSegmentSource->GetOutput();

	vtkSmartPointer<vtkCentrelinePartitioner> centrelinePartitioner = vtkSmartPointer<vtkCentrelinePartitioner>::New();
	centrelinePartitioner->SetInputData(vesselCentreline);
	centrelinePartitioner->SetPartitionLength(50);
	centrelinePartitioner->Update();

	vtkSmartPointer<vtkPolyData> centrelineSegments = centrelinePartitioner->GetOutput();
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputData(centrelineSegments);
	writer->SetFileName("TEST_227_ResampledBifurcation_Partitioned.vtk");
	//writer->SetFileName("TEST_721_ResampledBifurcation_Partitioned.vtk");
	writer->SetFileTypeToASCII();
	writer->Write();
	centrelinePartitioner->Print(std::cout);



	return EXIT_SUCCESS;
}

