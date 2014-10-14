#include <stdlib.h>
#include <assert.h>

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

#include <vtkAxesActor.h>
#include <vtkTransform.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>
#include <vtkPointData.h>

#include <vtkCallbackCommand.h>
#include <vtkCommand.h>

#include "vtkDbiharPatchFilter.h"

void PrintPoint(double *point)
{
	std::cout << "["<< point[0] << "," << point[1] << "," << point[2] << "]";
}

void KeypressCallbackFunction(vtkObject* caller, long unsigned int vtkNotUsed(eventId), void *clientData, void* vtkNotUsed(callData))
{
	vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);

	std::cout << "Key pressed: " << iren->GetKeySym() << std::endl;

	vtkSmartPointer<vtkActor> derivativesActor = reinterpret_cast<vtkActor*>(clientData);
	if(derivativesActor != 0 && std::string(iren->GetKeySym()) == "v")
	{
		derivativesActor->SetVisibility(!derivativesActor->GetVisibility());
		iren->Render();
	}
}

void showPolyData(vtkPolyData *input, vtkStructuredGrid *output)
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

		vtkDataArray *derivatives = input->GetPointData()->GetVectors(vtkDbiharPatchFilter::DERIV_ARR_NAME);
		if(derivatives != 0)
		{
			derivativesPresent = true;

			input->GetPointData()->SetActiveVectors(vtkDbiharPatchFilter::DERIV_ARR_NAME);

			vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
			arrowSource->Update();

			vtkSmartPointer<vtkGlyph3D> derivativesFilter = vtkSmartPointer<vtkGlyph3D>::New();
			derivativesFilter->SetInputData(input);
			derivativesFilter->SetScaleFactor(0.1);
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
	renderer->SetBackground(0.3, 0.6, 0.3); // Background color green

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
