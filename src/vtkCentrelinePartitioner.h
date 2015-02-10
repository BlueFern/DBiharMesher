#ifndef __vtkCentrelinePartitioner_h_
#define __vtkCentrelinePartitioner_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkPriorityQueue.h>

class vtkIdList;

/**
 * This filter partitions given vtkPolyData into segments using a user set partition
 * length as a guide. This bound should represent how many points (roughly) are to be
 * included to or from a bifurcation point to avoid having non-smooth boundaries between
 * branches. Therefore the segment size for a spine over such a bifurcation will be
 * closer to twice the partition length.
 *
 * The size of the straight segments between bifurcations is less important. Their
 * size will be roughly the input bound, and always greater than the minimum number
 * of points the dbihar patch filter requires to build an edge (so long as the input
 * data has enough points in each cell).
 *
 * The size of the segments is only ever rough due to the requirement for each segment
 * to have odd length. Using the given bound, the program attempts to divide a cell/branch
 * into a number of segments of that size (odd length). As a result of this process
 * the sizes of spines not participating in bifurcations tend to be less than or equal to
 * the partition length.
 *
 * \param vtkPolyData As is centreline data.
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
