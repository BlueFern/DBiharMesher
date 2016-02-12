#ifndef __vtkDbiharPatchSmooth_h_
#define __vtkDbiharPatchSmooth_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>

/**
 * vtkDbiharPatchSmooth is a filter, generated to smooth the arrangement of quads of an simplyfied
 * coronary artery bifurcation model (such as quadMeshFullc4080.vtp). To do so the algorithm must be
 * supplied with the original model as input. It then generates two new daugther branches, appends
 * them to the existing trunk and gives that new model as output. It creates new daughter branches
 * by taking their contours and applying the DbiharPatchFilter on each of them.
 * To be able to calculate the necessary data from the original model, the number of radial
 * quads ('numRadialQuads', halved) must be set to the same value as it is in the pipeline that
 * creates the input data (i.e. meshGen4080.cxx). With the aim of this parameter the algorithm is
 * able to calculate all the necessary data by itself. First it gets important coordinates (point in
 * the centre of bifurcation and endpoints of daughter branches) and calculates the model's properties
 * (lengths, radii, angles). Knowing that there are four endothelial cells per quad in axial direction
 * and the lengths of the branches, the number of axial quads per branch ('numAxialQuads0',
 * 'numAxialQuads1', 'numAxialQuads2') can be determined. In order to define the contour of a single
 * daugther branch, consecutive points are appended into an array. To get the desired coordinates specific
 * functions (several for-loops) are applied on the original vtkPolyData. After the points are set, they
 * are equipped with derivatives and the 'DbiharPatchFilter' implemented. As soon as both daughter
 * branches are created and attached to the parent branch, the model is been given as output.
 */


class vtkIdList;

class vtkDbiharPatchSmooth : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkDbiharPatchSmooth,vtkPolyDataAlgorithm);
	void PrintSelf(ostream& os, vtkIndent indent);

	static vtkDbiharPatchSmooth *New();

	vtkSetMacro(NumRadialQuads, int);

protected:
	vtkDbiharPatchSmooth();
	~vtkDbiharPatchSmooth() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);


private:
	vtkDbiharPatchSmooth(const vtkDbiharPatchSmooth&); // Not implemented.
	void operator=(const vtkDbiharPatchSmooth&); // Not implemented.

	int NumRadialQuads;


};

#endif
