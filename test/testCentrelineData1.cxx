
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

#include "wrapDbiharConfig.h"
#include "vtkCentrelineData.h"
#include "showPolyData.h"

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

// Get line ids.
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

void GetLineDirection(vtkPolyData *polyData, vtkIdType lineId, double *vector, direction_loc location)
{
	std::cout << __FUNCTION__ << std::endl;

	vtkSmartPointer<vtkIdList> line = GetLineIds(polyData, lineId);

	double p0[3];
	double p1[3];
	// TODO: Use error macro to check there are enough points for this.
	assert(line->GetNumberOfIds() > 0);
	if(location == d_start)
	{
		std::cout << "Getting (S) point with id " << line->GetId(1) << std::endl;
		polyData->GetPoint(line->GetId(0), p0);
		polyData->GetPoint(line->GetId(1), p1);
	}
	else if(location == d_end)
	{
		std::cout << "Getting (E) point with id " << line->GetId(line->GetNumberOfIds() - 2) << std::endl;
		polyData->GetPoint(line->GetId(line->GetNumberOfIds() - 2), p0);
		polyData->GetPoint(line->GetId(line->GetNumberOfIds() - 1), p1);
	}
	vtkMath::Subtract(p1, p0, vector);
}

void GetLineDirection(vtkPolyData *polyData, vtkIdType lineId, vtkIdType pointId, double *vector)
{
	vtkSmartPointer<vtkIdList> line = GetLineIds(polyData, lineId);

	// TODO: Use error macro to check there are enough points for this.
	assert(pointId - 1 >= 0);
	assert(line->GetNumberOfIds() > pointId);

	double p0[3];
	double p1[3];

	// Mapping from a local id in the line to polydata global point id.
	polyData->GetPoint(line->GetId(pointId - 1), p0);
	polyData->GetPoint(line->GetId(pointId), p1);

	vtkMath::Subtract(p1, p0, vector);
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
	//double tmp[3];
	vtkMath::Cross(c0, v0, c1);
	vtkMath::Cross(v1, c1, c1);
	// vtkMath::Normalize(c1);
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

	vtkSmartPointer<vtkCellArray> lines = resampledVesselCentreline->GetLines();
	lines->InitTraversal();

	// Obtain bifurcations from the centreline tree.
	std::map<vtkIdType, std::vector<vtkIdType> > bifurcations;
	// Remember the vessel segments not ending in bifurcations.
	vtkSmartPointer<vtkIdList> terminals = vtkSmartPointer<vtkIdList>::New();

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

		if(neighbourCellIds->GetNumberOfIds() > 0)
		{
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
		else
		{
			terminals->InsertNextId(lineId);
		}
	}

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

	// For first bifurcation.
	// Find cross product 0.
	// Translate it to start of step.
	// Find cross product 1, store it.
	// Translate cross product to start of step, label it 0.
	// Find cross product 1, store it.

	// 0 is for upward radii.
	vtkSmartPointer<vtkDoubleArray> radiiVectors0 = vtkSmartPointer<vtkDoubleArray>::New();
	radiiVectors0->SetNumberOfComponents(3);
	radiiVectors0->SetNumberOfTuples(resampledVesselCentreline->GetNumberOfPoints());
	radiiVectors0->SetName("radii_vectors");
	for(int i = 0; i < radiiVectors0->GetNumberOfTuples(); i++)
	{
		radiiVectors0->SetTuple3(i, 0, 0, 0);
	}

	double zero[3] = {0};
	double v0[3];
	double v1[3];
	//double v2[3];
	double c0[3];
	double c1[3];
	//double c2[3];

	// Keeping track of the average vectors at bifurcations.
	// This is a somewhat clunky way of doing so, which a map from
	// bifurcation point ids to indices in the bifurcation vectors array.
	std::map<vtkIdType, vtkIdType> lineIdToAverageVectorMap;
	vtkSmartPointer<vtkDoubleArray> bifurcationVectors = vtkSmartPointer<vtkDoubleArray>::New();
	bifurcationVectors->SetNumberOfComponents(3);

	for(std::map<vtkIdType, std::vector<vtkIdType> >::iterator it = bifurcations.begin(); it != bifurcations.end(); ++it)
	{
		GetLineDirection(resampledVesselCentreline, it->second[0], v0, d_start);
		GetLineDirection(resampledVesselCentreline, it->second[1], v1, d_start);

		// c0 orthogonal to v0 and v1.
		vtkMath::Cross(v0, v1, c0);
		vtkMath::Normalize(c0);

		// Global point id.
		vtkIdType bifId = GetPointIdFromLine(resampledVesselCentreline, it->first, d_end);

		// Set the length of c0 to the radius value at this point.
		vtkMath::MultiplyScalar(c0, resampledVesselCentreline->GetPointData()->GetScalars()->GetTuple1(bifId));

		// Store c0.
		radiiVectors0->SetTuple(bifId, c0);

		// Process the rest of the line.
		vtkSmartPointer<vtkIdList> lineIds = GetLineIds(resampledVesselCentreline, it->first);
		vtkMath::Normalize(c0);

		// v0 here is an average of the two vectors.
		vtkMath::Add(v0, v1, v0);
		vtkMath::MultiplyScalar(v0, 0.5);

		// Remember the average vector at the start of this bifurcation.
		lineIdToAverageVectorMap[it->second[0]] = bifurcationVectors->InsertNextTuple(v0);
		lineIdToAverageVectorMap[it->second[1]] = bifurcationVectors->InsertNextTuple(v0);

		// Rotation angle for radius at every point, to make sure there is a
		// gradual twist between start and end radius vector for each segment.
		double twistAngle = 0.0;
		double twistDirection = 1;

		// Figure out the twist angle between bifurcations.
		if(it->first > bifurcations.begin()->first)
		{
			vtkIdType parentBifId = GetPointIdFromLine(resampledVesselCentreline, it->first, d_start);
			bifurcationVectors->GetTuple((lineIdToAverageVectorMap[it->first]), v1);

			DoubleCross(v0, c0, v1, c1);
			vtkMath::Normalize(c1);

			double tmp[3];
			// c0 from parent bifurcation.
			radiiVectors0->GetTuple(parentBifId, tmp);

			// Figure out the angle between this bifurcation, and the closest upstearm bifurcation.
			twistAngle = vtkMath::DegreesFromRadians(AngleBetweenVectors(c1, tmp));
			std::cout << twistAngle << std::endl;

			// Figure out twist direction.
			vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
			transform->RotateWXYZ(twistAngle, v1);
			transform->TransformPoint(tmp, tmp);

			if(vtkMath::DegreesFromRadians(AngleBetweenVectors(tmp, c1)) < 1.0e-3)
			{
				twistDirection = -1;
			}
		}

		// Insert radii for the rest of the line.
		for(vtkIdType pId = lineIds->GetNumberOfIds() - 2; pId >= 1; pId--)
		{
			GetLineDirection(resampledVesselCentreline, it->first, pId, v1);

			DoubleCross(v0, c0, v1, c1);
			vtkMath::Normalize(c1);

			vtkIdType globalPointId = GetPointIdFromLine(resampledVesselCentreline, it->first, pId);
			vtkMath::MultiplyScalar(c1, resampledVesselCentreline->GetPointData()->GetScalars()->GetTuple1(globalPointId));

			// Add the twist.
			if(it->first > bifurcations.begin()->first)
			{
				double angleInc = twistAngle / (double)pId;
				vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
				transform->RotateWXYZ(twistDirection * angleInc, v1);
				transform->TransformPoint(c1, c1);
				twistAngle -= angleInc;
			}

			// Store at point pId vector data.
			radiiVectors0->SetTuple(globalPointId, c1);

			// Save v1 as v0 for next iteration.
			vtkMath::Add(v1, zero, v0);

			// Save c1 as c0 for next iteration.
			vtkMath::Add(c1, zero, c0);
		}
		// TODO: Need to insert first point for root trunk.
	}

	for(std::map<vtkIdType, vtkIdType>::iterator it = lineIdToAverageVectorMap.begin(); it != lineIdToAverageVectorMap.end(); ++it)
	{
		std::cout << it->first << " => " << it->second << std::endl;
		double tuple[3];
		bifurcationVectors->GetTuple(it->second, tuple);
		PrintPoint(tuple); std::cout << std::endl;
	}



	for(vtkIdType id = 0; id < terminals->GetNumberOfIds(); id++)
	{
		vtkIdType lId = terminals->GetId(id);

		// Process the rest of the line.
		vtkSmartPointer<vtkIdList> lineIds = GetLineIds(resampledVesselCentreline, lId);

		vtkIdType parentBifId = GetPointIdFromLine(resampledVesselCentreline, lId, d_start);
		std::cout << parentBifId << std::endl;

		GetLineDirection(resampledVesselCentreline, lId, 1, v0);
		NegateVector(v0, v0);
		// c0 from parent bifurcation.
		radiiVectors0->GetTuple(parentBifId, c0);

		for(vtkIdType pId = 2; pId < lineIds->GetNumberOfIds() - 1; pId++)
		{
			GetLineDirection(resampledVesselCentreline, lId, pId + 1, v1);
			NegateVector(v1, v1);

			DoubleCross(v0, c0, v1, c1);
			vtkMath::Normalize(c1);

			vtkIdType globalPointId = GetPointIdFromLine(resampledVesselCentreline, lId, pId);
			vtkMath::MultiplyScalar(c1, resampledVesselCentreline->GetPointData()->GetScalars()->GetTuple1(globalPointId));

			// Store at point pId vector data.
			radiiVectors0->SetTuple(globalPointId, c1);

			// Save v1 as v0 for next iteration.
			vtkMath::Add(v1, zero, v0);

			// Save c1 as c0 for next iteration.
			vtkMath::Add(c1, zero, c0);
		}
	}

	resampledVesselCentreline->GetPointData()->SetVectors(radiiVectors0);
	resampledVesselCentreline->GetPointData()->SetActiveVectors("radii_vectors");

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

#if 1
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

