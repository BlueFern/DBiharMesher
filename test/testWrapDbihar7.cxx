#include <sstream>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkIdList.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkAppendPoints.h>

#include "vtkDbiharStatic.h"
#include "vtkDbiharPatchFilter.h"
#include "wrapDbiharConfig.h"

#define SSTR( x ) dynamic_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()

int main(int argc, char* argv[]) {

	std::cout << "Starting " << __FILE__ << std::endl;

	vtkSmartPointer<vtkAppendPoints> appendFilter = vtkSmartPointer<vtkAppendPoints>::New();

	for(int i = 0; i < 3; i++)
	{
		vtkSmartPointer<vtkXMLPolyDataReader> patchReader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
		patchReader->SetFileName(SSTR(std::string(TEST_DATA_DIR) + "/inputPatch_" << i << ".vtp").c_str());
		patchReader->Update();

		vtkSmartPointer<vtkDbiharPatchFilter> patchFilter = vtkSmartPointer<vtkDbiharPatchFilter>::New();

		// Set the bounds of the UV space.
		patchFilter->SetA(0.0);
		patchFilter->SetB(0.7);
		patchFilter->SetC(0.0);
		patchFilter->SetD(1.2);

		// Set the dimensions conditions.
		patchFilter->SetMQuads(28);
		patchFilter->SetNQuads(88);

		// Set solution method.
		patchFilter->SetIFlag(4);

		patchFilter->SetInputData(patchReader->GetOutput());
		patchFilter->Update();

		appendFilter->AddInputData(patchFilter->GetOutput());
	}

	appendFilter->Update();

	vtkSmartPointer<vtkPolyData> outputPoints = appendFilter->GetOutput();
	vtkSmartPointer<vtkIdList> vertices = vtkSmartPointer<vtkIdList>::New();
	for(int i = 0; i < outputPoints->GetNumberOfPoints(); i++)
	{
		vertices->InsertNextId(i);
	}
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	cells->InsertNextCell(vertices);
	outputPoints->SetVerts(cells);

#if 0
	vtkSmartPointer<vtkXMLPolyDataWriter> outputWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	outputWriter->SetFileName("outputPatchPointsCombined.vtp");
	outputWriter->SetInputData(appendFilter->GetOutput());
	outputWriter->Update();
#endif

	vtkDbiharStatic::ShowPolyData(outputPoints);

	std::cout << "Exiting " << __FILE__ << std::endl;

	return EXIT_SUCCESS;
}

