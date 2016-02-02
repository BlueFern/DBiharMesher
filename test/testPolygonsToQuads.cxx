#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkSubdivideQuadFilter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkIdList.h>
#include <vtkCell.h>
#include <vtkGeometryFilter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>

#include "wrapDbiharConfig.h"
#include "vtkPolygonsToQuads.h"
#include "vtkDbiharStatic.h"



int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	reader->SetFileName((std::string(TEST_DATA_DIR) + "/test/ec_data_t_25.vtu").c_str());
	reader->Update();

	vtkSmartPointer<vtkGeometryFilter> geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
	geometryFilter->SetInputData(reader->GetOutput());
	geometryFilter->Update();

	vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
	pd->DeepCopy(geometryFilter->GetOutput());


	vtkSmartPointer<vtkPolygonsToQuads> polygonsToQuads = vtkSmartPointer<vtkPolygonsToQuads>::New();
	polygonsToQuads->SetInputData(pd);
	polygonsToQuads->Update();

	vtkDbiharStatic::ShowPolyData(polygonsToQuads->GetOutput());
#if 1
	vtkDbiharStatic::WritePolyData(polygonsToQuads->GetOutput(), "polygonsToQuadsTest.vtp");
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}
