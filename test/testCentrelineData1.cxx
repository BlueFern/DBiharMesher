
#include <map>
#include <list>

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
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkMath.h>
#include <vtkTransform.h>

#include <vtkXMLPolyDataWriter.h>
#include <vtkGenericDataObjectWriter.h>

#include "vtkCentrelineData.h"
#include "showPolyData.h"

#include "wrapDbiharConfig.h"

double returnPoint[3] = {0.0};
void NegateVector(const double v[3], double nV[3])
{
	nV[0] = -v[0];
	nV[1] = -v[1];
	nV[2] = -v[2];
}
double* NegateVector(const double v[3])
{
	NegateVector(v, returnPoint);
	return returnPoint;
}

typedef enum {d_start = 0, d_end = 1} direction_loc;

// Get line ids for a given line.
vtkSmartPointer<vtkIdList> GetLineIds(vtkPolyData *polyData, vtkIdType lineId)
{
	bool lineFound = false;
	vtkSmartPointer<vtkCellArray> lines = polyData->GetLines();
	vtkSmartPointer<vtkIdList> line = vtkSmartPointer<vtkIdList>::New();
	lines->InitTraversal();
	for(vtkIdType id = 0; id < lines->GetNumberOfCells(); id++)
	{
		lines->GetNextCell(line);
		if(id == lineId)
		{
			lineFound = true;
			break;
		}
	}
	if(lineFound == false)
	{
		// TODO: Use warning macro.
		std::cerr << "WARNING: Line " << lineId << " not found." << std::endl;
		line->SetNumberOfIds(0);
	}
	return line;
}

// Get global point id.
vtkIdType GetPointIdFromLine(vtkPolyData *polyData, vtkIdType lineId, direction_loc location)
{
	vtkSmartPointer<vtkIdList> line = GetLineIds(polyData, lineId);

	// TODO: Use warning macro to check there are enough points for this.
	if(line->GetNumberOfIds() == 0)
	{
		std::cerr << "WARNING: Line " << lineId << " is empty." << std::endl;
	}

	vtkIdType pId = -1;
	if(location == d_start)
	{
		pId = line->GetId(0);
	}
	else if(location == d_end)
	{
		pId = line->GetId(line->GetNumberOfIds() - 1);
	}
	return pId;
}

// Get global point id.
vtkIdType GetPointIdFromLine(vtkPolyData *polyData, vtkIdType lineId, vtkIdType pointId)
{
	vtkSmartPointer<vtkIdList> line = GetLineIds(polyData, lineId);

	// TODO: Use error macro to check there are enough points for this.
	assert(line->GetNumberOfIds() > pointId);

	vtkIdType pId = -1;

	// Mapping from a local id in the line to polydata global point id.
	pId = line->GetId(pointId);

	return pId;
}

// Get the vector indicating the direction of the centreline at this location.
void GetLineDirection(vtkPolyData *polyData, vtkIdType lineId, double *vector, direction_loc location)
{
	vtkSmartPointer<vtkIdList> line = GetLineIds(polyData, lineId);

	double p0[3];
	double p1[3];

	// TODO: Use VTK error macro to check there are enough points for this.
	assert(line->GetNumberOfIds() > 0);

	if(location == d_start)
	{
		//std::cout << "Getting (S) point with id " << line->GetId(1) << std::endl;
		polyData->GetPoint(line->GetId(1), p0);
		polyData->GetPoint(line->GetId(0), p1);
	}
	else if(location == d_end)
	{
		//std::cout << "Getting (E) point with id " << line->GetId(line->GetNumberOfIds() - 2) << std::endl;
		polyData->GetPoint(line->GetId(line->GetNumberOfIds() - 1), p0);
		polyData->GetPoint(line->GetId(line->GetNumberOfIds() - 2), p1);
	}
	vtkMath::Subtract(p1, p0, vector);
}

// Get the vector indicating the direction of the centreline at this location.
void GetLineDirection(vtkPolyData *polyData, vtkIdType lineId, vtkIdType pointId, double *vector)
{
	vtkSmartPointer<vtkIdList> line = GetLineIds(polyData, lineId);

	// TODO: Use error macro to check there are enough points for this.
	assert(pointId >= 0);
	assert(line->GetNumberOfIds() > pointId);

	double p0[3];
	double p1[3];

	// Mapping from a local id in the line to polydata global point id.
	if(pointId != 0)
	{
		polyData->GetPoint(line->GetId(pointId), p0);
		polyData->GetPoint(line->GetId(pointId - 1), p1);
		vtkMath::Subtract(p1, p0, vector);
		//GetLineDirection(polyData, line->GetId(pointId), vector);
	}
	else
	{
		// This is a bit fishy, because the direction at the fist point is calculated an the direction at the next point.
		GetLineDirection(polyData, lineId, pointId + 1, vector);
	}
}

// This is taken from vtkMath class in VTK Nightly in October 2014.
double AngleBetweenVectors(const double v1[3], const double v2[3])
{
  double cross[3];
  vtkMath::Cross(v1, v2, cross);
  return atan2(vtkMath::Norm(cross), vtkMath::Dot(v1, v2));
}

void DoubleCross(const double v0[3], const double c0[3], const double v1[3], double c1[3])
{
	vtkMath::Cross(c0, v0, c1);
	vtkMath::Cross(v1, c1, c1);
}

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

	vtkSmartPointer<vtkIdList> verts = vtkSmartPointer<vtkIdList>::New();
	for(vtkIdType vId = 0; vId < resampledVesselCentreline->GetNumberOfPoints(); vId++)
	{
		verts->InsertNextId(vId);
	}

	vtkSmartPointer<vtkCellArray> vertsArray = vtkSmartPointer<vtkCellArray>::New();
	vertsArray->InsertNextCell(verts);

	resampledVesselCentreline->SetVerts(vertsArray);

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

#if 0
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(resampledVesselCentreline);
	//mapper->SetInputData(vesselCentreline);

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(512, 512);
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkInteractorStyleSwitch *iss = vtkInteractorStyleSwitch::SafeDownCast(renderWindowInteractor->GetInteractorStyle());
	iss->SetCurrentStyleToTrackballCamera();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	renderer->AddActor(actor);

	renderWindow->Render();
	renderWindowInteractor->Start();
	renderWindow->Finalize();
#endif

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}

