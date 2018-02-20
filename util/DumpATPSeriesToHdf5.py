"""
Reads in a number of .vtp ATP files and writes them to hdf5 files. Each time-step
exists as a dataset within the hdf5 file for the particular branch of the input
geometry.
"""
import h5py
import os
import vtk
import glob

originalTimeStep = 0.01
timeStep = 0.01
taskMeshIn = ''
atpMeshPattern = ''
axialQuads = 0
circQuads = 0
numECsPerCol = 4
numSMCsPerRow = 4

numECsPerRow = numSMCsPerRow * 5
numSMCsPerCol = numECsPerCol * 13

numECsPerQuad = numECsPerRow * numECsPerCol
numSMCsPerQuad = numSMCsPerCol * numSMCsPerRow

atpHdf5Files = [
"files/parent_atp.h5",
"files/left_daughter_atp.h5",
"files/right_daughter_atp.h5",
]

def writeHdf5():
    
    if timeStep < 0.01:
        exit("Timestep is too small, choose 0.01 or larger")
    
    # This is where the data is for testing purposes.
    print("Current working directory:", os.getcwd())

    numQuadsPerRing = circQuads

    # Working with the task mesh junt to figure out the quads and rows numbers.
    taskMeshReader = vtk.vtkXMLPolyDataReader()
    taskMeshReader.SetFileName(taskMeshIn)
    taskMeshReader.Update()

    taskMesh = taskMeshReader.GetOutput()

    # Get the range of branch labels.
    labelRange = [0, 0]
    taskMesh.GetCellData().GetScalars().GetRange(labelRange, 0)

    # Convert label range to a list of labels.
    labelRange = range(int(labelRange[0]), int(labelRange[1]) + 1)
    print("Labels found in task mesh:", labelRange)

    # Store the number of rings for each label. 
    numRingsPerLabel = {}   

    # For every label in the range of labels we want to extract all cells/quads.
    for label in labelRange:
        
        # Use this filter to extract the cells for a given label value.
        branchSelector = vtk.vtkThreshold()
        branchSelector.SetInputData(taskMesh)
        branchSelector.ThresholdBetween(label,label);
        branchSelector.Update()

        taskMeshBranch = branchSelector.GetOutput()

        numQuadRowsPerBranch = taskMeshBranch.GetNumberOfCells() / numQuadsPerRing;
        numRingsPerLabel[label] = numQuadRowsPerBranch
        
    atpFiles = sorted(glob.glob(atpMeshPattern))

    parentFile = h5py.File(atpHdf5Files[0], 'w')
    leftBranchFile = h5py.File(atpHdf5Files[1], 'w')
    rightBranchFile = h5py.File(atpHdf5Files[2], 'w')
    
    for atpIndex in range(0, len(atpFiles), int(timeStep / originalTimeStep)):
        print("Time step" + str(atpIndex * timeStep))
        print("Reading", atpFiles[atpIndex], "at index", atpIndex)

        atpMeshReader = vtk.vtkXMLPolyDataReader()
        atpMeshReader.SetFileName(atpFiles[atpIndex])
        atpMeshReader.Update()
    
        atpMesh = atpMeshReader.GetOutput()

        # For every label in the range of labels we want to extract all ECs.
        for label in labelRange:
    
            # Keep track of how many branches we need to skip.
            numECsPerLabel = numQuadsPerRing * numRingsPerLabel[label] * numECsPerQuad
            atpCellOffset = label * numECsPerLabel
    
            # Collect cell ids to select.
            selectionIds = vtk.vtkIdTypeArray()
            for sId in range(0, int(numECsPerLabel)):
                selectionIds.InsertNextValue(int(atpCellOffset) + sId)
    
            # Create selecion node.
            selectionNode = vtk.vtkSelectionNode()
            selectionNode.SetFieldType(selectionNode.CELL)
            selectionNode.SetContentType(selectionNode.INDICES)
            selectionNode.SetSelectionList(selectionIds)
    
            # Create selection.
            selection = vtk.vtkSelection()
            selection.AddNode(selectionNode)
    
            # Use vtkSelection filter.
            selectionExtractor = vtk.vtkExtractSelection()
            selectionExtractor.SetInputData(0, atpMesh)
            selectionExtractor.SetInputData(1, selection)
            selectionExtractor.Update()
    
            extractedCells = selectionExtractor.GetOutput()

            # Ring ids list for traversal.
            ringIds = range(0, int(numRingsPerLabel[label]))
            ringIds = list(ringIds)
            ringIds.reverse()

            # Number of ECs rows is the number of ECs per quad.
            rowIds = range(0, numECsPerCol)
            rowIds = list(rowIds)
            rowIds.reverse()

            # Decide which h5 files to write to.
            pointsOf = ''
            
            if label == 0:
                pointsOf = parentFile
            elif label == 1:
                pointsOf = leftBranchFile
            elif label == 2:
                pointsOf = rightBranchFile

            dset = pointsOf.create_dataset("/" + str(atpIndex), (numECsPerLabel,), 'f')
            i = 0
    
            # Iterate over the rings in reverse order.
            for ringNum in ringIds:
                # Iterate over the 'imaginary' quads of cells in normal order.
                for quadNum in range(0, numQuadsPerRing):
                    # Iterate over the rows of cells in reverse order.
                    # Calculate the 'real' id for the 'imaginary' quad.
                    quadId = ringNum * numQuadsPerRing + quadNum
                    # Iterate over rows of cells in reverse order.
                    for rowNum in rowIds:
                        # Iterate over the rows of cells in normal order.
                        for cellNum in range(0, numECsPerRow):
                            # Calculate the 'real' ec cell id and get the corresponding cell.
                            realId = quadId * numECsPerQuad + rowNum * numECsPerRow + cellNum

                            atpVal = extractedCells.GetCellData().GetArray("ATP").GetValue(realId)

                            # Insert the value into the dataset.
                            dset[i] = atpVal
                            i += 1
    parentFile.close()
    leftBranchFile.close()
    rightBranchFile.close()

    print("All done ...")
    
def main():
    print("This script is to be run with global parameters (input, output files, etc.) set in the calling script.")

if __name__ == '__main__':
    print("Starting", os.path.basename(__file__))
    main()
    print("Exiting", os.path.basename(__file__))