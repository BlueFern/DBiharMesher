#include <stdlib.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkGenericDataObjectReader.h>

#include "wrapDbiharConfig.h"
#include "vtkScalarRadiiToVectorsFilter.h"
#include "vtkSkipSegmentFilter.h"
#include "vtkDbiharStatic.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/testSkipSegmentFilter1.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(vesselCentreline);
	scalarRadiiToVectorsFilter->Update();

	vtkSmartPointer<vtkSkipSegmentFilter> skipSegmentFilter = vtkSmartPointer<vtkSkipSegmentFilter>::New();
	skipSegmentFilter->SetInputData(scalarRadiiToVectorsFilter->GetOutput());
	skipSegmentFilter->SetInlet(true);
	//skipSegmentFilter->SetOutlet(true);
	skipSegmentFilter->SetSkipSize(16);
	skipSegmentFilter->SetPointId(0);
	//skipSegmentFilter->SetPointId(504);
	//skipSegmentFilter->SetPointId(658);
	skipSegmentFilter->SetNumberOfRadialQuads(28);
	skipSegmentFilter->Update();

	vtkSmartPointer<vtkPolyData> centrelineSegments = skipSegmentFilter->GetOutput();

	vtkDbiharStatic::WritePolyData(centrelineSegments, "skipSegmentTest.vtp");

	skipSegmentFilter->Print(std::cout);

	return EXIT_SUCCESS;
}
