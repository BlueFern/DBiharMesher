# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 17:09:57 2015

This script read an ECs mesh produced by WrapDbihar code and converts it to
a ATP profile map for the legacy Coupled Cells code.

The code reorders cells and produces the required files in TXT format.
VTK files are written out for visual verification.
"""
import h5py
import os
import vtk


taskMeshIn = ''
ecMeshIn = ''
atpMeshIn = ''
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

# EC VTP files used for their geometry in visual verification of ATP mesh.
ecVTPFiles = [
"vtk/ec_mesh_parent.vtp",
"vtk/ec_mesh_left_daughter.vtp",
"vtk/ec_mesh_right_daughter.vtp",
]

def writeHdf5():
    # This is where the data is for testing purposes.
    print("Current working directory:", os.getcwd())
    print(taskMeshIn)
    
    numQuadsPerRing = circQuads

    # Working with the task mesh junt to figure out the quads and rows numbers.
    taskMeshReader = vtk.vtkXMLPolyDataReader()
    taskMeshReader.SetFileName(taskMeshIn)
    taskMeshReader.Update()

    taskMesh = taskMeshReader.GetOutput()
    print(taskMesh.GetNumberOfPoints())
    
    ecMeshReader = vtk.vtkXMLPolyDataReader()
    ecMeshReader.SetFileName(ecMeshIn)
    ecMeshReader.Update()
    
    ecMesh = ecMeshReader.GetOutput()
    print(ecMesh.GetNumberOfPoints())
    
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

    # Working with EC mesh only
    atpMeshReader = vtk.vtkXMLPolyDataReader()
    atpMeshReader.SetFileName(atpMeshIn)
    atpMeshReader.Update()

    # Original ECs mesh to work with.
    atpMesh = atpMeshReader.GetOutput()
    print("There are", atpMesh.GetNumberOfCells(), "ATP values in total ...")

    parentFile = h5py.File(atpHdf5Files[0], 'w')
    leftBranchFile = h5py.File(atpHdf5Files[1], 'w')
    rightBranchFile = h5py.File(atpHdf5Files[2], 'w')
    
    appendPolyData = vtk.vtkAppendPolyData();

    # For every label in the range of labels we want to extract all ECs.
    for label in labelRange:
        
        ecMeshReader = vtk.vtkXMLPolyDataReader()
        ecMeshReader.SetFileName(ecVTPFiles[label])
        ecMeshReader.Update()
        tmpPolyData = ecMeshReader.GetOutput()

        # Keep track of how many branches we need to skip.
        numECsPerLabel = numQuadsPerRing * numRingsPerLabel[label] * numECsPerQuad
        atpCellOffset = label * numECsPerLabel

        print("atpCellOffset", atpCellOffset)

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
        
        reorderedATPArray = vtk.vtkDoubleArray()
        reorderedATPArray.SetName("initialATP")

        # Decide which TXT files to write to.
        pointsOf = ''
        
        if label == 0:
            pointsOf = parentFile
        elif label == 1:
            pointsOf = leftBranchFile
        elif label == 2:
            pointsOf = rightBranchFile

        print("Writing H5 file for ECs ATP:")
        print(pointsOf)
        dset = pointsOf.create_dataset("/atp", (numECsPerLabel,), 'f')
        
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
                        
                        atpVal = extractedCells.GetCellData().GetArray("initialATP").GetValue(realId)
                        
                        reorderedATPArray.InsertNextValue(atpVal)
                        

                        # Write the value to the txt file.
                        dset[i] = atpVal
                        i += 1
        
        tmpPolyData.GetCellData().SetScalars(reorderedATPArray)
        appendPolyData.AddInputData(tmpPolyData)
                        

    parentFile.close()
    leftBranchFile.close()
    rightBranchFile.close()
    
    print("Writing reorderd ATP map for verification...")
    appendPolyData.Update()
    reorderedATPWriter = vtk.vtkXMLPolyDataWriter()
    reorderedATPWriter.SetInputData(appendPolyData.GetOutput())
    reorderedATPWriter.SetFileName("vtk/reordered_atp.vtp")
    reorderedATPWriter.Update()
    

    print("All done ...")
    
def main():
    print("This script is to be run with global parameters (input, output files, etc.) set in the calling script.")

if __name__ == '__main__':
    print("Starting", os.path.basename(__file__))
    main()
    print("Exiting", os.path.basename(__file__))
