#include <vtkSmartPointer.h>
#include <vtkObjectFactory.h>

#include "vtkCentrelineData.h"

vtkStandardNewMacro(vtkCentrelineData);

vtkCentrelineData::vtkCentrelineData()
{
	std::cout << __FUNCTION__ << std::endl;;
}

int vtkCentrelineData::GetNumVessels()
{
	return 0;
}

int vtkCentrelineData::GetNumBifurcations()
{
	return 0;
}

void vtkCentrelineData::PrintSelf(ostream &os, vtkIndent indent)
{
	this->Superclass::PrintSelf(os, indent);
}
