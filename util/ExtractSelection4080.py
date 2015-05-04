import os
import re
import glob
import vtk

def sortNicely(l): 
  """ Sort the given list in the way that humans expect. 
  """ 
  convert = lambda text: int(text) if text.isdigit() else text 
  alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
  l.sort( key=alphanum_key )


# Prepare selection id array.
start = 6240
numRingsPerBranch = 34
numQuadsPerRing = 40
numSMCsPerCol = 4 * 13

ringOffset = numQuadsPerRing * numSMCsPerCol * 4
branchOffset = numRingsPerBranch * ringOffset

selectionIds = vtk.vtkIdTypeArray()

for branch in [0, 2]:
	thisBranchOffset = branch * branchOffset

	for ring in range(34):
		thisRingOffset = ring * ringOffset

		for cell in range(numSMCsPerCol):
			thisOffset = start + thisBranchOffset + thisRingOffset + cell
			selectionIds.InsertNextValue(thisOffset)

# Create selection node.
selectionNode = vtk.vtkSelectionNode()
selectionNode.SetFieldType(selectionNode.CELL)
selectionNode.SetContentType(selectionNode.INDICES)
selectionNode.SetSelectionList(selectionIds)

# Create selection.
selection = vtk.vtkSelection()
selection.AddNode(selectionNode)

# Get the list of files to process and sort it.
fileList = glob.glob('solutionVTU/smc_Data_t_*.vtu')
sortNicely(fileList)

# Process every file by extracting selection.
for inFile in fileList:

	print 'Reading file', inFile

	reader = vtk.vtkXMLUnstructuredGridReader()
	reader.SetFileName(inFile)
	
	print '\tExtracting ', selectionIds.GetNumberOfTuples(), 'cells...'
	selectionExtractor = vtk.vtkExtractSelection()
	selectionExtractor.SetInputConnection(reader.GetOutputPort())
	selectionExtractor.SetInputData(1, selection)

	baseName = os.path.basename(inFile)
	number = [int(s) for s in re.split('([0-9]+)', baseName) if s.isdigit()]
	outFile = 'selection/' + 'extracted_SMC_Line_' + str(number[0]) + '.vtu'

	print '\t\tSaving file', outFile

	writer = vtk.vtkXMLUnstructuredGridWriter()
	writer.SetInputConnection(selectionExtractor.GetOutputPort())
	writer.SetFileName(outFile)
	writer.Update()

	print '\n'

print 'All done...'

