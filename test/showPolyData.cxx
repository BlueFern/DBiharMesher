#include <stdlib.h>
#include <assert.h>

#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkAxesActor.h>
#include <vtkTransform.h>

void showPolyData(vtkPolyData *input, vtkPolyData *output)
{
	// At least one should be not zero.
	assert(input != 0 || output != 0);

	vtkSmartPointer<vtkActor> inputActor;
	if(input != NULL)
	{
		vtkSmartPointer<vtkPolyDataMapper> inputMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		inputMapper->SetInputData(input);

		inputActor = vtkSmartPointer<vtkActor>::New();
		inputActor->SetMapper(inputMapper);
		inputActor->GetProperty()->SetColor(0,1,0);
	}

	vtkSmartPointer<vtkActor> outputActor;
	if(output != 0)
	{
		vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
		for(int pId = 0; pId < output->GetNumberOfPoints(); pId++)
		{
			idList->InsertNextId(pId);
		}

		vtkSmartPointer<vtkCellArray> verts = vtkSmartPointer<vtkCellArray>::New();
		verts->InsertNextCell(idList);
		output->SetVerts(verts);

		vtkSmartPointer<vtkPolyDataMapper> outputMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		outputMapper->SetInputData(output);

		outputActor = vtkSmartPointer<vtkActor>::New();
		outputActor->SetMapper(outputMapper);
		outputActor->GetProperty()->SetPointSize(2);
		outputActor->GetProperty()->SetColor(1,0,0);
	}

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkInteractorStyleSwitch *iss = vtkInteractorStyleSwitch::SafeDownCast(renderWindowInteractor->GetInteractorStyle());
	iss->SetCurrentStyleToTrackballCamera();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	if(inputActor != 0)
	{

		renderer->AddActor(inputActor);
	}

	if(outputActor != 0)
	{
		renderer->AddActor(outputActor);
	}

	vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->Translate(-15.0, -20.0, 0.0);
	axes->SetUserTransform(transform);
	renderer->AddActor(axes);

	renderWindow->Render();
	renderWindowInteractor->Start();

	renderWindow->Finalize();
}
