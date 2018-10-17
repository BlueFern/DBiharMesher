#include <cmath>
#include <algorithm>
#include <sstream>

#include <time.h>

#include <vtkMath.h>
#include <vtkIdList.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkSTLWriter.h>
#include <vtkPNGWriter.h>

#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkProperty.h>
#include <vtkPolyDataMapper.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSetMapper.h>

#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>

#include <vtkWindowToImageFilter.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>

#include <vtkCallbackCommand.h>
#include <vtkCommand.h>

#include "vtkDbiharStatic.h"

const char *vtkDbiharStatic::RADII_VECTORS_ARR_NAME = {"radiiVectors"};
const char *vtkDbiharStatic::RADII_SCALARS_ARR_NAME = {"radiiScalars"};
const char *vtkDbiharStatic::DERIV_ARR_NAME = {"derivVectors"};
const char *vtkDbiharStatic::BRANCH_ID_ARR_NAME = {"branchId"};
const char *vtkDbiharStatic::GRID_COORDS_ARR_NAME = {"gridCoords"};

// The following constants are in micrometres.
const double vtkDbiharStatic::SMC_CIRC = 50 / 2.0;
const double vtkDbiharStatic::SMC_AXIAL = 5 / 2.0;
const double vtkDbiharStatic::EC_CIRC = 10 / 2.0;
const double vtkDbiharStatic::EC_AXIAL = 65 / 2.0;

void KeypressCallbackFunction(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void *clientData, void* vtkNotUsed(callData))
{
	vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);

	std::cout << "Key pressed: " << iren->GetKeySym() << std::endl;

	vtkSmartPointer<vtkActor> actor = reinterpret_cast<vtkActor*>(clientData);
	if(actor != 0 && std::string(iren->GetKeySym()) == "v")
	{
		actor->SetVisibility(!actor->GetVisibility());
		iren->Render();
	}

	if(std::string(iren->GetKeySym()) == "s")
	{
		vtkSmartPointer<vtkWindowToImageFilter> windowToImageFilter = vtkSmartPointer<vtkWindowToImageFilter>::New();
		windowToImageFilter->SetInput(iren->GetRenderWindow());
		windowToImageFilter->SetMagnification(3);
		windowToImageFilter->ReadFrontBufferOff(); // Read from the back buffer.
		windowToImageFilter->Update();

		time_t _tm =time(NULL);
		struct tm *curtime = localtime(&_tm);
		std::ostringstream outputFileNameStream;
		outputFileNameStream << "Screenshot_" << asctime(curtime) << ".png";

		std::string outputFileName = outputFileNameStream.str();
		outputFileName.erase(std::remove(outputFileName.begin(), outputFileName.end(), '\n'), outputFileName.end());
		std::replace(outputFileName.begin(), outputFileName.end(), ' ', '_');
		std::replace(outputFileName.begin(), outputFileName.end(), ':', '-');

		vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
		writer->SetFileName(outputFileName.c_str());
		writer->SetInputConnection(windowToImageFilter->GetOutputPort());
		writer->Write();
	}
}

std::string vtkDbiharStatic::GetTimeStamp()
{
	time_t rawtime;
	struct tm * timeinfo;
	char buffer[80];

	time (&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer,80,"%d-%m-%Y %I:%M:%S",timeinfo);
	std::string str(buffer);

	return str;
}

void vtkDbiharStatic::WritePolyData(vtkPolyData *input, std::string fileName)
{
	std::clog << vtkDbiharStatic::GetTimeStamp() << ": Writing file " << fileName << std::endl;
	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetInputData(input);
	writer->SetFileName(fileName.c_str());
	writer->Update();

}

void vtkDbiharStatic::WriteStlData(vtkPolyData *input, std::string fileName)
{
	std::clog << vtkDbiharStatic::GetTimeStamp() << ": Writing file " << fileName << std::endl;
	vtkSmartPointer<vtkSTLWriter> writer = vtkSmartPointer<vtkSTLWriter>::New();
	writer->SetInputData(input);
	writer->SetFileName(fileName.c_str());
	writer->Update();
}

void vtkDbiharStatic::ShowPolyDataWithGrid(vtkPolyData *input, vtkStructuredGrid *output, double derivateScaling)
{
	// At least one should be not zero.
	assert(input != 0 || output != 0);

	vtkSmartPointer<vtkActor> inputActor;
	vtkSmartPointer<vtkActor> derivativesActor;

	bool derivativesPresent = false;

	if(input != NULL)
	{
		vtkSmartPointer<vtkPolyDataMapper> inputMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		inputMapper->SetInputData(input);

		inputActor = vtkSmartPointer<vtkActor>::New();
		inputActor->SetMapper(inputMapper);
		inputActor->GetProperty()->SetColor(1,0,0);
		inputActor->GetProperty()->SetLineWidth(1.5);
		inputActor->GetProperty()->SetRepresentationToWireframe();

		vtkDataArray *derivatives = input->GetPointData()->GetVectors("derivVectors");
		if(derivatives != 0)
		{
			derivativesPresent = true;

			input->GetPointData()->SetActiveVectors("derivVectors");

			vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
			arrowSource->Update();

			vtkSmartPointer<vtkGlyph3D> derivativesFilter = vtkSmartPointer<vtkGlyph3D>::New();
			derivativesFilter->SetInputData(input);
			derivativesFilter->SetScaleFactor(derivateScaling);
			derivativesFilter->SetSourceConnection(arrowSource->GetOutputPort());
			derivativesFilter->SetScaleModeToScaleByVector();
			derivativesFilter->SetVectorModeToUseVector();
			derivativesFilter->OrientOn();
			derivativesFilter->Update();

			vtkSmartPointer<vtkPolyDataMapper> derivativesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
			derivativesMapper->SetInputConnection(derivativesFilter->GetOutputPort());

			derivativesActor = vtkSmartPointer<vtkActor>::New();
			derivativesActor->SetMapper(derivativesMapper);
		}
	}

	vtkSmartPointer<vtkActor> gridActor;
	if(output != 0)
	{
		vtkSmartPointer<vtkDataSetMapper> gridMapper = vtkSmartPointer<vtkDataSetMapper>::New();
		gridMapper->SetInputData(output);

		gridActor = vtkSmartPointer<vtkActor>::New();
		gridActor->SetMapper(gridMapper);
		gridActor->GetProperty()->EdgeVisibilityOn();
		gridActor->GetProperty()->SetEdgeColor(0, 0, 1);
		gridActor->GetProperty()->SetOpacity(0.7);
	}

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	// renderer->SetBackground(0.3, 0.6, 0.3); // Green background.
	renderer->SetBackground(0.317, 0.341, 0.431); // ParaView background.

	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(600, 600);
	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkInteractorStyleSwitch *iss = vtkInteractorStyleSwitch::SafeDownCast(renderWindowInteractor->GetInteractorStyle());
	iss->SetCurrentStyleToTrackballCamera();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	if(inputActor != 0)
	{

		renderer->AddActor(inputActor);
	}

	if(gridActor != 0)
	{
		renderer->AddActor(gridActor);
	}

	vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
	if(input != 0)
	{
		double bounds[6];
		input->GetBounds(bounds);
		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		// Translate the axis to lower left corner.
		transform->Translate(bounds[0] - (bounds[1] - bounds[0]) / 8.0, bounds[2] - (bounds[3] - bounds[2]) / 8.0, 0.0);
		axes->SetUserTransform(transform);
	}
	renderer->AddActor(axes);

	if(derivativesPresent)
	{
		renderer->AddActor(derivativesActor);

		vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		keypressCallback->SetCallback(KeypressCallbackFunction);
		keypressCallback->SetClientData(derivativesActor);
		renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);
	}

	renderWindow->Render();
	renderWindowInteractor->Start();

	renderWindow->Finalize();
}

void vtkDbiharStatic::ShowPolyData(vtkPolyData *input, double vectorScaling)
{
	vtkSmartPointer<vtkPolyDataMapper> inputMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	inputMapper->SetInputData(input);

	vtkSmartPointer<vtkActor> inputActor = vtkSmartPointer<vtkActor>::New();
	inputActor->SetMapper(inputMapper);
	//inputActor->GetProperty()->SetColor(1,0,0);
	inputActor->GetProperty()->SetLineWidth(1.5);
	inputActor->GetProperty()->SetPointSize(5);
	inputActor->GetProperty()->SetRepresentationToSurface();
	inputActor->GetProperty()->SetEdgeVisibility(true);

	vtkSmartPointer<vtkActor> derivativesActor = vtkSmartPointer<vtkActor>::New();
	vtkDataArray *derivatives = input->GetPointData()->GetVectors();
	if(derivatives != 0)
	{
		vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
		arrowSource->Update();

		vtkSmartPointer<vtkGlyph3D> derivativesFilter = vtkSmartPointer<vtkGlyph3D>::New();
		derivativesFilter->SetInputData(input);
		derivativesFilter->SetScaleFactor(vectorScaling);
		derivativesFilter->SetSourceConnection(arrowSource->GetOutputPort());
		derivativesFilter->SetScaleModeToScaleByVector();
		derivativesFilter->SetVectorModeToUseVector();
		derivativesFilter->OrientOn();
		derivativesFilter->Update();

		vtkSmartPointer<vtkPolyDataMapper> derivativesMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		derivativesMapper->SetInputConnection(derivativesFilter->GetOutputPort());

		derivativesActor->SetMapper(derivativesMapper);
	}

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	// renderer->SetBackground(0.3, 0.6, 0.3); // Green background.
	renderer->SetBackground(0.317, 0.341, 0.431); // ParaView background.

	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(600, 600);
	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkInteractorStyleSwitch *iss = vtkInteractorStyleSwitch::SafeDownCast(renderWindowInteractor->GetInteractorStyle());
	iss->SetCurrentStyleToTrackballCamera();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(inputActor);
	renderer->AddActor(derivativesActor);

	vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
	if(input != 0)
	{
		double bounds[6];
		input->GetBounds(bounds);
		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		// Translate the axis to lower left corner.
		transform->Translate(bounds[0] - (bounds[1] - bounds[0]) / 8.0, bounds[2] - (bounds[3] - bounds[2]) / 8.0, 0.0);
		axes->SetUserTransform(transform);
	}
	renderer->AddActor(axes);

	vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
	keypressCallback->SetCallback(KeypressCallbackFunction);
	keypressCallback->SetClientData(derivativesActor);
	renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);

	renderWindow->Render();
	renderWindowInteractor->Start();

	renderWindow->Finalize();
}

void vtkDbiharStatic::DoubleCross(const double v0[3], const double c0[3], const double v1[3], double c1[3])
{
	vtkMath::Cross(c0, v0, c1);
	vtkMath::Cross(v1, c1, c1);
}

void vtkDbiharStatic::PrintDataArray(vtkDataArray *dataArray)
{
	dataArray->Print(std::cout);
	for(int tupleNum = 0; tupleNum < dataArray->GetNumberOfTuples(); tupleNum++)
	{
		for(int componentNum = 0; componentNum < dataArray->GetNumberOfComponents(); componentNum++)
		{
			std::cout << dataArray->GetVariantValue(tupleNum * dataArray->GetNumberOfComponents() + componentNum).ToString() << " ";
		}
		std::cout << std::endl;
	}
}
