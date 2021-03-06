
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedIntArray.h>
#include <vtkCellArray.h>

#include <vtkPointData.h>
#include <vtkXMLStructuredGridReader.h>

#include "wrapDbiharConfig.h"
#include "vtkDbiharStatic.h"
#include "vtkCentrelinePartitioner.h"

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkGenericDataObjectReader> vesselCentrelineReader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
	vesselCentrelineReader->SetFileName((std::string(TEST_DATA_DIR) + "/test/testCentrelinePartitioner1_0.vtk").c_str());
	vesselCentrelineReader->Update();

	vtkPolyData *vesselCentreline = vtkPolyData::SafeDownCast(vesselCentrelineReader->GetOutput());


	vtkSmartPointer<vtkCentrelinePartitioner> centrelinePartitioner = vtkSmartPointer<vtkCentrelinePartitioner>::New();
	centrelinePartitioner->SetInputData(vesselCentreline);
	centrelinePartitioner->SetPartitionLength(50);

	vtkSmartPointer<vtkIdList> EndPoints = vtkSmartPointer<vtkIdList>::New();
	int x = 0;
	switch(x)
	{
		case 0:
			EndPoints->InsertNextId(220);
			EndPoints->InsertNextId(420); // End points on same branch.

			break;

		case 1:
			EndPoints->InsertNextId(60); // Starting from a different cell, will ignore one side of the tree.
			break;

		case 2:
			EndPoints->InsertNextId(30);
			EndPoints->InsertNextId(54);
			EndPoints->InsertNextId(90);
			break;

		case 3:
			EndPoints->InsertNextId(60);
			EndPoints->InsertNextId(1306); // Will incur a warning, 1306 will never be reached.
			break;

		case 4:
			EndPoints->InsertNextId(10);
			EndPoints->InsertNextId(60);
			EndPoints->InsertNextId(1370); // End points lie next to bifurcations.
			break;

		case 5:
			EndPoints->InsertNextId(1260); // One cell and only one specified endpoint.
			break;
	}


	//centrelinePartitioner->SetEndPoints(EndPoints);
	centrelinePartitioner->Update();


	centrelinePartitioner->Print(std::cout);

	vtkDbiharStatic::WritePolyData(centrelinePartitioner->GetOutput(), "partitionedCentreline.vtp");

	return EXIT_SUCCESS;
}

