#ifndef __vtkCentrelinePartitioner_h_
#define __vtkCentrelinePartitioner_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPriorityQueue.h>

class vtkIdList;

/**
 * This filter partitions given vtkPolyData into segments using a user set partition
 * length as a guide. The sizes of each segment created will roughly be the given
 * length. In rare cases, segments have the potential to be closer to half the length,
 * due to rounding.
 *
 * The size of the segments is only ever rough due to the requirement for each segment
 * to have odd length. Using the given bound, the program attempts to divide a cell/branch
 * into a number of segments of that size (odd length). As a result of this process
 * the actual length of each segment not participating in a bifurcation is described by:
 * (partition length / 2) < actual length <= partition length.
 *
 * Before partitioning bifurcation points are easily found by checking the last id in each cell,
 * but this information is lost when the partitioner segments the tree. Therefore the bifurcation
 * points are recorded and set as vertices in the output vtkPolyData structure (using vtkPolyVertex).
 *
 * For similar reasons end points are also recorded and stored in a vtkPolyVertex structure.
 *
 * \param vtkPolyData As centreline data.
 *
 * \param EndPoints An optional Id list that specifies points to build segments between.
 * The first point is where the partitioner should start, and all others are used as new
 * end points for the associated cell. The user is responsible for giving sensible end
 * points so that the cell sizes are still all odd. If points are given that belong
 * in cells that branched out before the given starting point they will be ignored.
 * If it is not set the entire centreline will be partitioned.
 *
 * \param PartitionLength A number that determines the rough size of the segments between bifurcations,
 * and the length to be traveled before and after a bifurcation in such segments. If the
 * number of points around a bifurcation is too small, the mesh is less likely to be smooth.
 *
 * \return vtkPolyData which has cells/lines as the partitioned segments, and points and
 * point data identical to the input.
 *
 */

class vtkCentrelinePartitioner : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkCentrelinePartitioner,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkCentrelinePartitioner *New();
	vtkSetMacro(PartitionLength, int);

	virtual void SetEndPoints(vtkIdList* EndPoints);

protected:
	vtkCentrelinePartitioner();
	~vtkCentrelinePartitioner() {};

	static void ProgressFunction(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData);

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkCentrelinePartitioner(const vtkCentrelinePartitioner&); // Not implemented.
	void operator=(const vtkCentrelinePartitioner&); // Not implemented.

	void joinIdLists(vtkSmartPointer<vtkIdList> previous, vtkSmartPointer<vtkIdList> current,
					   vtkSmartPointer<vtkIdList> joined);
	void reverseIdList(vtkSmartPointer<vtkIdList> spine, vtkSmartPointer<vtkIdList> reversedSpine);

	int PartitionLength;
	vtkSmartPointer<vtkIdList> EndPoints;
};

#endif
