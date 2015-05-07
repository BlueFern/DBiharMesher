import os
import re
import vtk

def sortNicely(l): 
    """ Sort the given list in the way that humans expect. 
    """ 
    convert = lambda text: int(text) if text.isdigit() else text 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    l.sort( key=alphanum_key )

# These parameters are to be set in the calling script.
start = 0
numRingsPerBranch = 0
numQuadsPerRing = 0
outputPattern = ''
branches = []
numSMCsPerCol = 4 * 13

def ExtractSelection(fileList):
    # Report our CWD just for testing purposes.
    print "CWD:", os.getcwd()

    ringOffset = numQuadsPerRing * numSMCsPerCol * 4
    branchOffset = numRingsPerBranch * ringOffset

    # Prepare selection id array.    
    selectionIds = vtk.vtkIdTypeArray()
    
    for branch in branches:
        thisBranchOffset = branch * branchOffset
        # print thisBranchOffset,
    
        for ring in range(numRingsPerBranch):
            thisRingOffset = ring * ringOffset
            # print thisRingOffset, 
            
            for cell in range(numSMCsPerCol):
                thisOffset = start + thisBranchOffset + thisRingOffset + cell
                # print thisOffset,
                selectionIds.InsertNextValue(thisOffset)
    
    # Create selection node.
    selectionNode = vtk.vtkSelectionNode()
    selectionNode.SetFieldType(selectionNode.CELL)
    selectionNode.SetContentType(selectionNode.INDICES)
    selectionNode.SetSelectionList(selectionIds)
    
    # Create selection.
    selection = vtk.vtkSelection()
    selection.AddNode(selectionNode)
    
    sortNicely(fileList)
    
    # Process every file by extracting selection.
    for inFile in fileList:
        print 'Reading file', inFile
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(inFile)
        
        print '\tExtracting ', selectionIds.GetNumberOfTuples(), 'cells...'
        selectionExtractor = vtk.vtkExtractSelection()
        selectionExtractor.SetInputConnection(reader.GetOutputPort())
        if vtk.vtkVersion().GetVTKMajorVersion() > 5:
            selectionExtractor.SetInputData(1, selection)
        else:
            selectionExtractor.SetInput(1, selection)
        
        baseName = os.path.basename(inFile)
        number = [int(s) for s in re.split('([0-9]+)', baseName) if s.isdigit()]
        outFile = outputPattern + str(number[0]) + '.vtu'
        
        print '\t\tSaving file', outFile
        writer = vtk.vtkXMLUnstructuredGridWriter()
        writer.SetInputConnection(selectionExtractor.GetOutputPort())
        writer.SetFileName(outFile)
        writer.Update()
        
        print '\n'
    
    print 'All done...'

def Usage():
	print "This script is to be run with global parameters (mesh dimensions, output file name pattern, etc.) set in the calling script."
	
if __name__ == '__main__':
	print "Starting", os.path.basename(__file__)
	Usage()
	print "Exiting", os.path.basename(__file__)

