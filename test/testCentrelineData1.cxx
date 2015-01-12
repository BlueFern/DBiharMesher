
#include <map>
#include <list>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectReader.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>

#include <vtkIdList.h>
#include <vtkCellArray.h>

#include <vtkXMLPolyDataWriter.h>
#include <vtkGenericDataObjectWriter.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>

#include "vtkCentrelineData.h"
#include "showPolyData.h"

#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/227A_Centreline.vtk").c_str());
	//vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/721A_Centreline.vtk").c_str());
	//vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/SyntheticBifurcation_1.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	vtkSmartPointer<vtkCentrelineData> centrelineSegmentSource = vtkSmartPointer<vtkCentrelineData>::New();
	centrelineSegmentSource->DebugOn();
	centrelineSegmentSource->SetEdgeLength(4 * 65e-3);
	centrelineSegmentSource->SetCentrelineData(vesselCentreline);


	std::cout << 4 * 65e-3 << ", " << 65e-3 * 4 << std::endl;

	vtkPolyData *resampledVesselCentreline = centrelineSegmentSource->GetOutput();

	std::cout << "Number of input points: " << vesselCentreline->GetNumberOfPoints() << std::endl;
	std::cout << "Number of input lines: " << vesselCentreline->GetNumberOfLines() << std::endl;

	std::cout << "Number of output points: " << resampledVesselCentreline->GetNumberOfPoints() << std::endl;
	std::cout << "Number of output lines: " << resampledVesselCentreline->GetNumberOfLines() << std::endl;

	vtkSmartPointer<vtkIdList> verts = vtkSmartPointer<vtkIdList>::New();
	for(vtkIdType vId = 0; vId < resampledVesselCentreline->GetNumberOfPoints(); vId++)
	{
		verts->InsertNextId(vId);
	}

	vtkSmartPointer<vtkCellArray> vertsArray = vtkSmartPointer<vtkCellArray>::New();
	vertsArray->InsertNextCell(verts);

	resampledVesselCentreline->SetVerts(vertsArray);

	showPolyData1(resampledVesselCentreline, 0.5);

#if 1
	vtkSmartPointer<vtkXMLPolyDataWriter> tmpWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	tmpWriter->SetInputData(resampledVesselCentreline);
	tmpWriter->SetFileName("227A_CentrelineResampled_4ECs.vtp");
	//tmpWriter->SetFileName("721A_CentrelineResampled_4ECs.vtp");
	//tmpWriter->SetFileName("resampledSyntheticBifurcation_1.vtp");
	tmpWriter->Write();

	vtkSmartPointer<vtkGenericDataObjectWriter> writer = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetInputData(resampledVesselCentreline);
	writer->SetFileName("227A_CentrelineResampled_4ECs.vtk");
	//writer->SetFileName("721A_CentrelineResampled_4ECs.vtk");
	//writer->SetFileName("resampledSyntheticBifurcation_1.vtk");
	writer->Write();
#endif

#if 0
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(resampledVesselCentreline);
	//mapper->SetInputData(vesselCentreline);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(512, 512);
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkInteractorStyleSwitch *iss = vtkInteractorStyleSwitch::SafeDownCast(renderWindowInteractor->GetInteractorStyle());
	iss->SetCurrentStyleToTrackballCamera();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderer->AddActor(actor);

	renderWindow->Render();
	renderWindowInteractor->Start();
	renderWindow->Finalize();
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}

