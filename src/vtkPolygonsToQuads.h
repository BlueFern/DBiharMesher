#ifndef __vtkPolygonsToQuads_h_
#define __vtkPolygonsToQuads_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

/*
 * This filter accepts meshes and converts all 6-sided polygons to quads. Cell data is
 * copied (both quads created from any given polygon have identical cell data).
 *
 */
class vtkPolygonsToQuads : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkPolygonsToQuads,vtkPolyDataAlgorithm);
	static vtkPolygonsToQuads *New();

protected:
	vtkPolygonsToQuads();
	~vtkPolygonsToQuads() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);
	void copyPointsArray(double array1[], double array2[]);

private:
	vtkPolygonsToQuads(const vtkPolygonsToQuads&); // Not implemented.
	void operator=(const vtkPolygonsToQuads&); // Not implemented.

};

#endif
