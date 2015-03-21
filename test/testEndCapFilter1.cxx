#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>

#include "vtkEndCapFilter.h"
#include "wrapDbiharConfig.h"
#include "vtkDbiharStatic.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
	reader->SetFileName((std::string(TEST_DATA_DIR) + "/endCapInputData_1.vtp").c_str());
	reader->Update();

	vtkSmartPointer<vtkEndCapFilter> endCapFilter = vtkSmartPointer<vtkEndCapFilter>::New();
	endCapFilter->SetInputData(reader->GetOutput());
	endCapFilter->Update();

	vtkDbiharStatic::ShowPolyDataWithGrid(endCapFilter->GetOutput(), NULL);

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
