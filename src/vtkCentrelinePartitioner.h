/*
 * vtkCentrelinePartioner.h
 *
 *  Created on: 18/12/2014
 *      Author: sed59
 */

#ifndef __vtkCentrelinePartitioner_h_
#define __vtkCentrelinePartitioner_h_

#include <vtkAlgorithm.h>
#include <vtkPolyDataAlgorithm.h>

class vtkIdList;

class vtkCentrelinePartitioner : public vtkPolyDataAlgorithm {
public:
	vtkTypeMacro(vtkCentrelinePartitioner,vtkPolyDataAlgorithm);

	static vtkCentrelinePartitioner *New();

	vtkSetMacro(Bound, int);

protected:
	vtkCentrelinePartitioner();
	~vtkCentrelinePartitioner() {};

	int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);

private:
	vtkCentrelinePartitioner(const vtkCentrelinePartitioner&); // Not implemented.
	void operator=(const vtkCentrelinePartitioner&); // Not implemented.

	int Bound;

};




#endif
