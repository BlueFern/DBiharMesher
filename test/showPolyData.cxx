#include <stdlib.h>
#include <assert.h>
#include <time.h>

#include <algorithm>
#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSetMapper.h>
#include <vtkProperty.h>

#include <vtkIdList.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkPointData.h>

#include <vtkCallbackCommand.h>
#include <vtkCommand.h>

#include "showPolyData.h"

void PrintPoint(double *point)
{
	std::cout << "["<< point[0] << "," << point[1] << "," << point[2] << "]";
}

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

void showPolyData(vtkPolyData *input, vtkStructuredGrid *output, double derivateScaling)
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
	renderer->SetBackground(0.3, 0.6, 0.3); // Green background.

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

void showPolyData1(vtkPolyData *input, double vectorScaling)
{
	vtkSmartPointer<vtkPolyDataMapper> inputMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	inputMapper->SetInputData(input);

	vtkSmartPointer<vtkActor> inputActor = vtkSmartPointer<vtkActor>::New();
	inputActor->SetMapper(inputMapper);
	inputActor->GetProperty()->SetColor(1,0,0);
	inputActor->GetProperty()->SetLineWidth(1.5);
	inputActor->GetProperty()->SetPointSize(5);

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
	renderer->SetBackground(0.3, 0.6, 0.3); // Green background.

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

void showGrids(std::vector<vtkSmartPointer<vtkStructuredGrid> > grids, vtkSmartPointer<vtkPolyData> centreline)
{

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(0.3, 0.6, 0.3); // Green background.

	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(600, 600);
	renderWindow->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkInteractorStyleSwitch *iss = vtkInteractorStyleSwitch::SafeDownCast(renderWindowInteractor->GetInteractorStyle());
	iss->SetCurrentStyleToTrackballCamera();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	if(centreline != NULL)
	{
		vtkSmartPointer<vtkPolyDataMapper> centrelineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		centrelineMapper->SetInputData(centreline);

		vtkSmartPointer<vtkActor> inputActor = vtkSmartPointer<vtkActor>::New();
		inputActor->SetMapper(centrelineMapper);
		inputActor->GetProperty()->SetColor(1,0,0);
		inputActor->GetProperty()->SetLineWidth(1.5);
		inputActor->GetProperty()->SetPointSize(5);

		vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
		arrowSource->Update();

		vtkSmartPointer<vtkGlyph3D> arrowsFilter = vtkSmartPointer<vtkGlyph3D>::New();
		arrowsFilter->SetInputData(centreline);
		arrowsFilter->SetScaleFactor(1.0);
		arrowsFilter->SetSourceConnection(arrowSource->GetOutputPort());
		arrowsFilter->SetScaleModeToScaleByVector();
		arrowsFilter->SetVectorModeToUseVector();
		arrowsFilter->OrientOn();
		arrowsFilter->Update();

		vtkSmartPointer<vtkPolyDataMapper> arrowsMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		arrowsMapper->SetInputConnection(arrowsFilter->GetOutputPort());

		vtkSmartPointer<vtkActor> arrowsActor = vtkSmartPointer<vtkActor>::New();
		arrowsActor->SetMapper(arrowsMapper);

		renderer->AddActor(inputActor);
		renderer->AddActor(arrowsActor);

		double bounds[6];
		centreline->GetBounds(bounds);
		vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
		// Translate the axis to lower left corner.
		transform->Translate(bounds[0] - (bounds[1] - bounds[0]) / 8.0, bounds[2] - (bounds[3] - bounds[2]) / 8.0, 0.0);
		vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
		axes->SetUserTransform(transform);
		renderer->AddActor(axes);

		vtkSmartPointer<vtkCallbackCommand> keypressCallback = vtkSmartPointer<vtkCallbackCommand>::New();
		keypressCallback->SetCallback(KeypressCallbackFunction);
		keypressCallback->SetClientData(arrowsActor);
		renderWindowInteractor->AddObserver(vtkCommand::KeyPressEvent, keypressCallback);
	}

	for(int i = 0; i < grids.size(); i++)
	{
		vtkSmartPointer<vtkDataSetMapper> gridMapper = vtkSmartPointer<vtkDataSetMapper>::New();
		gridMapper->SetInputData(grids[i]);
		vtkSmartPointer<vtkActor> gridActor = vtkSmartPointer<vtkActor>::New();
		gridActor->SetMapper(gridMapper);
		gridActor->GetProperty()->EdgeVisibilityOn();
		gridActor->GetProperty()->SetEdgeColor(0, 0, 1);
		gridActor->GetProperty()->SetOpacity(1.0);
		renderer->AddActor(gridActor);
	}

	renderWindow->Render();
	renderWindowInteractor->Start();

	renderWindow->Finalize();
}
