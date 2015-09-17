
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
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>

#include "vtkCentrelineResampler.h"
#include "vtkDbiharStatic.h"

#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/testCentrelineResampler0.vtk").c_str());
	//vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/testCentrelineResampler1.vtk").c_str());
	//vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/testCentrelineResamplerBifurcation.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	vtkSmartPointer<vtkCentrelineResampler> centrelineSegmentSource = vtkSmartPointer<vtkCentrelineResampler>::New();
	centrelineSegmentSource->DebugOn();
	centrelineSegmentSource->SetEdgeLength(4 * 65e-3);
	centrelineSegmentSource->SetInputData(vesselCentreline);


	std::cout << 4 * 65e-3 << ", " << 65e-3 * 4 << std::endl;
	centrelineSegmentSource->Update();
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

	vtkDbiharStatic::ShowPolyData(resampledVesselCentreline, 0.5);
	centrelineSegmentSource->Print(std::cout);

#if 1
	vtkDbiharStatic::WritePolyData(resampledVesselCentreline, "227A_CentrelineResampled_4ECs.vtp");
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

