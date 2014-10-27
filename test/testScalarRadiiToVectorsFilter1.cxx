
#include <map>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectReader.h>

#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleSwitch.h>

#include <vtkCellArray.h>
#include <vtkPointData.h>

#include <vtkXMLPolyDataWriter.h>
#include <vtkGenericDataObjectWriter.h>

#include "vtkCentrelineData.h"
#include "vtkScalarRadiiToVectorsFilter.h"
#include "showPolyData.h"

#include "wrapDbiharConfig.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/227A_Centreline.vtk").c_str());
	//vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/721A_Centreline.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());

	vtkSmartPointer<vtkCentrelineData> centrelineSegmentSource = vtkSmartPointer<vtkCentrelineData>::New();
	centrelineSegmentSource->SetCentrelineData(vesselCentreline);

	vtkPolyData *resampledVesselCentreline = centrelineSegmentSource->GetOutput();

	std::cout << "Number of input points: " << vesselCentreline->GetNumberOfPoints() << std::endl;
	std::cout << "Number of input lines: " << vesselCentreline->GetNumberOfLines() << std::endl;

	std::cout << "Number of output points: " << resampledVesselCentreline->GetNumberOfPoints() << std::endl;
	std::cout << "Number of output lines: " << resampledVesselCentreline->GetNumberOfLines() << std::endl;

	vtkSmartPointer<vtkCellArray> lines = resampledVesselCentreline->GetLines();
	lines->InitTraversal();

	// Obtain bifurcations from the centreline tree.
	std::map<vtkIdType, std::vector<vtkIdType> > bifurcations; // TODO: This needs a more appropriate name.

	// Find bifurcations and terminals.
	for(vtkIdType lineId = 0; lineId < lines->GetNumberOfCells(); lineId++)
	{
		vtkSmartPointer<vtkIdList> lineIds = vtkSmartPointer<vtkIdList>::New();
		lines->GetNextCell(lineIds);
		vtkSmartPointer<vtkIdList> lastId = vtkSmartPointer<vtkIdList>::New();
		lastId->InsertNextId(lineIds->GetId(lineIds->GetNumberOfIds() - 1));

		// WARNING: The code in the GetCellNeighbors method changes the traversal location in the vtkCellArray. In my view it is utterly retarded.
		vtkIdType traverseLocation = resampledVesselCentreline->GetLines()->GetTraversalLocation();

		vtkSmartPointer<vtkIdList> neighbourCellIds = vtkSmartPointer<vtkIdList>::New();
		resampledVesselCentreline->GetCellNeighbors(lineId, lastId, neighbourCellIds);

		// Restore traversal location.
		resampledVesselCentreline->GetLines()->SetTraversalLocation(traverseLocation);

		std::cout << "Line " << lineId << " shares last point with " << neighbourCellIds->GetNumberOfIds() << " cells." << std::endl;

		std::vector<vtkIdType> idList;

		// The assumption is that non-bifurcating segments are not divided in segments.
		if(neighbourCellIds->GetNumberOfIds() == 1)
		{
			// vtkWarningMacro("Found a point shared between only two cells, which breaks the assumption of non-bifurcating segments integrity.");
		}

		for(vtkIdType pId = 0; pId < neighbourCellIds->GetNumberOfIds(); pId++)
		{
			idList.push_back(neighbourCellIds->GetId(pId));
		}
		bifurcations[lineId] = idList;
	}

	// Print the map.
	for(std::map<vtkIdType, std::vector<vtkIdType> >::iterator it = bifurcations.begin(); it != bifurcations.end(); ++it)
	{
		std::vector<vtkIdType> v = it->second;
		std::cout << it->first << " => ";
		for(std::vector<vtkIdType>::iterator i = v.begin(); i != v.end(); ++i)
		{
			std::cout << *i << " ";
		}
		std::cout << endl;
	}

	vtkSmartPointer<vtkScalarRadiiToVectorsFilter> scalarRadiiToVectorsFilter = vtkSmartPointer<vtkScalarRadiiToVectorsFilter>::New();
	scalarRadiiToVectorsFilter->SetInputData(resampledVesselCentreline);
	scalarRadiiToVectorsFilter->Update();

	showPolyData1(resampledVesselCentreline, 0.5);

#if 1
	vtkSmartPointer<vtkXMLPolyDataWriter> tmpWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	tmpWriter->SetInputData(resampledVesselCentreline);
	tmpWriter->SetFileName("resampledCentreline.vtp");
	tmpWriter->Write();

	vtkSmartPointer<vtkGenericDataObjectWriter> writer = vtkSmartPointer<vtkGenericDataObjectWriter>::New();
	writer->SetFileName("resampledCentreline.vtk");
	writer->SetInputData(resampledVesselCentreline);
	writer->Write();
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}

