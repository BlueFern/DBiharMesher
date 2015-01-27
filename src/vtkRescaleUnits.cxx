/**
 * Rescales a data set (as vtkPolyData) by a specified amount. This includes both the points
 * and point data.
 */

#include <algorithm>

#include <vtkObjectFactory.h>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkDataObject.h>
#include <vtkSmartPointer.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>

#include "vtkRescaleUnits.h"

vtkStandardNewMacro(vtkRescaleUnits);

vtkRescaleUnits::vtkRescaleUnits()
{
	this->SetNumberOfInputPorts(1);
	this->SetNumberOfOutputPorts(1);

	this->Scale = 1;
}

int vtkRescaleUnits::RequestData(vtkInformation *vtkNotUsed(request),
		vtkInformationVector **inputVector, vtkInformationVector *outputVector)
{
	// Get the input and output.
	vtkPolyData* input = vtkPolyData::GetData(inputVector[0], 0);
	vtkPolyData* output = vtkPolyData::GetData(outputVector, 0);

	vtkSmartPointer<vtkTransformFilter> transformFilter = vtkSmartPointer<vtkTransformFilter>::New();
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();

	vtkSmartPointer<vtkDoubleArray> doubleArray = vtkSmartPointer<vtkDoubleArray>::New();
	double scale[3] = {this->Scale, this->Scale, this->Scale}; // Used as diagonal values in scalar matrix.
	transform->Scale(scale);

	transformFilter->SetInputData(input);
	transformFilter->SetTransform(transform);
	transformFilter->Update(); // Point data is unchanged, however.

	doubleArray->SetNumberOfValues(input->GetPointData()->GetNumberOfTuples());
	doubleArray = vtkDoubleArray::SafeDownCast(input->GetPointData()->GetScalars());

	// Multiply each point data by given scale.
	for (int i = 0; i < doubleArray->GetNumberOfTuples(); i++)
	{
		doubleArray->SetValue(i, doubleArray->GetValue(i) * this->Scale);
	}

	output->GetPointData()->SetScalars(doubleArray);
	output->ShallowCopy(transformFilter->GetOutput());

	return 1;
}

void vtkRescaleUnits::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
	os << indent << "Given scale: " << this->Scale << "\n";
}
