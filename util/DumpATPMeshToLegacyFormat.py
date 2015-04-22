# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 17:09:57 2015

This script read an ECs mesh produced by WrapDbihar code and converts it to
a ATP profile map for the legacy Coupled Cells code.

The code reorders cells and produces the required files in TXT format.
VTK files are written out for visual verification.
"""

import os
import vtk

numECsPerCol = 4
numSMCsPerRow = 4

numECsPerRow = numSMCsPerRow * 5
numSMCsPerCol = numECsPerCol * 13

numECsPerQuad = numECsPerRow * numECsPerCol
numSMCsPerQuad = numSMCsPerCol * numSMCsPerRow

'''
numQuadsPerRing0 = 12
taskMeshIn = "quadMeshFullc216.vtp"
ecMeshIn = "quadMeshFullECc216.vtp"
atpMeshIn = "quadMeshFullATPc216.vtp"
''' and None

'''
numQuadsPerRing0 = 48
taskMeshIn = "quadMeshFullc4032.vtp"
ecMeshIn = "quadMeshFullECc4032.vtp"
atpMeshIn = "quadMeshFullATPc4032.vtp"
''' and None

'''
numQuadsPerRing0 = 40
taskMeshIn = "quadMeshFullc4080.vtp"
ecMeshIn = "quadMeshFullECc4080.vtp"
atpMeshIn = "quadMeshFullATPc4080.vtp"
''' and None

# '''
numQuadsPerRing0 = 64
taskMeshIn = "quadMeshFullc8064.vtp"
ecMeshIn = "quadMeshFullECc8064.vtp"
atpMeshIn = "quadMeshFullATPc8064.vtp"
# ''' and None


# VTK files to write.
atpVTKFiles = [
"vtk/atp_parent.vtp",
"vtk/atp_left_daughter.vtp",
"vtk/atp_right_daughter.vtp",
]

atpTXTFiles = [
"files/parent_atp.txt",
"files/left_daughter_atp.txt",
"files/right_daughter_atp.txt",
]

def main():
    # This is where the data is for testing purposes.
    print "Current working directory:", os.getcwd()

    # Working with the task mesh junt to figure out the quads and rows numbers.
    taskMeshReader = vtk.vtkXMLPolyDataReader()
    taskMeshReader.SetFileName(taskMeshIn)
    taskMeshReader.Update()

    taskMesh = taskMeshReader.GetOutput()
    print taskMesh.GetNumberOfPoints()
    
    ecMeshReader = vtk.vtkXMLPolyDataReader()
    ecMeshReader.SetFileName(ecMeshIn)
    ecMeshReader.Update()
    
    ecMesh = ecMeshReader.GetOutput()
    print ecMesh.GetNumberOfPoints()
    
    # Get the range of branch labels.
    labelRange = [0, 0]
    taskMesh.GetCellData().GetScalars().GetRange(labelRange, 0)

    # Convert label range to a list of labels.
    labelRange = range(int(labelRange[0]), int(labelRange[1]) + 1)
    print "Labels found in task mesh:", labelRange

    # Store the number of rings for each label. 
    numRingsPerLabel = {}   

    # For every label in the range of labels we want to extract all cells/quads.
    for label in labelRange:
        
        # Use this filter to extract the cells for a given label value.
        branchSelector = vtk.vtkThreshold()
        branchSelector.SetInput(taskMesh)
        branchSelector.ThresholdBetween(label,label);
        branchSelector.Update()

        taskMeshBranch = branchSelector.GetOutput()

        numQuadRowsPerBranch = taskMeshBranch.GetNumberOfCells() / numQuadsPerRing0;
        numRingsPerLabel[label] = numQuadRowsPerBranch

    # Working with EC mesh only
    atpMeshReader = vtk.vtkXMLPolyDataReader()
    atpMeshReader.SetFileName(atpMeshIn)
    atpMeshReader.Update()

    # Original ECs mesh to work with.
    atpMesh = atpMeshReader.GetOutput()
    print "There are", atpMesh.GetNumberOfCells(), "ATP values in total ..."

    # Prepare the stupid TXT files for output.    
    parentFile = open(atpTXTFiles[0], 'w')
    leftBranchFile = open(atpTXTFiles[1], 'w')
    rightBranchFile = open(atpTXTFiles[2], 'w')

    # For every label in the range of labels we want to extract all ECs.
    for label in labelRange:

        # Keep track of how many branches we need to skip.
        numECsPerLabel = numQuadsPerRing0 * numRingsPerLabel[label] * numECsPerQuad
        atpCellOffset = label * numECsPerLabel

        print "atpCellOffset", atpCellOffset

        # Collect cell ids to select.
        selectionIds = vtk.vtkIdTypeArray()
        for sId in range(0, numECsPerLabel):
            selectionIds.InsertNextValue(atpCellOffset + sId)

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
        selectionExtractor.SetInput(0, atpMesh)
        selectionExtractor.SetInput(1, selection)
        selectionExtractor.Update()

        extractedCells = selectionExtractor.GetOutput()

        # Ring ids list for traversal.
        ringIds = range(0, numRingsPerLabel[label])
        ringIds.reverse()

        # Number of ECs rows is the number of ECs per quad.
        rowIds = range(0, numECsPerCol)
        rowIds.reverse()

        # New vtkPoints for storing reordered points.
        reorderedPoints = vtk.vtkPoints()

        # New vtkCellArray for storing reordeced cells.
        reorderedCellArray = vtk.vtkCellArray()
        
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

        print "Writing TXT file for ECs ATP:"
        print pointsOf

        # Iterate over the rings in reverse order.
        for ringNum in ringIds:
            # Iterate over the 'imaginary' quads of cells in normal order.
            for quadNum in range(0, numQuadsPerRing0):
                # Iterate over the rows of cells in reverse order.
                # Calculate the 'real' id for the 'imaginary' quad.
                quadId = ringNum * numQuadsPerRing0 + quadNum
                # Iterate over rows of cells in reverse order.
                for rowNum in rowIds:
                    # Iterate over the rows of cells in normal order.
                    for cellNum in range(0, numECsPerRow):
                        # Calculate the 'real' ec cell id and get the corresponding cell.
                        realId = quadId * numECsPerQuad + rowNum * numECsPerRow + cellNum
                        # Get the corresponding cell from the ecMesh.
                        cell = ecMesh.GetCell(label * numECsPerLabel + realId)
                        cellPoints = cell.GetPoints()
                        
                        # This is for writing out the surface files for visualisation.
                        newCell = vtk.vtkQuad()
                        for ptId in range(0, cellPoints.GetNumberOfPoints()):
                            newCell.GetPointIds().SetId(ptId, reorderedPoints.InsertNextPoint(cellPoints.GetPoint(ptId)))
                        
                        reorderedCellArray.InsertNextCell(newCell)
                        
                        atpVal = extractedCells.GetCellData().GetArray("initialATP").GetValue(realId)
                        
                        reorderedATPArray.InsertNextValue(atpVal)

                        # Write the value to the txt file.
                        pointsOf.write(format(atpVal, '.6f') + '\n')
                        
        # Create new vtkPolyData object for the new reordered mesh.
        reorderedATPBranch = vtk.vtkPolyData()

        # Insert our new points.
        reorderedATPBranch.SetPoints(reorderedPoints)

        # Set the reordered cells to the reordered ATP mesh.
        reorderedATPBranch.SetPolys(reorderedCellArray)
        
        # Set the reordered ATP values to the reordered ATP mesh.
        reorderedATPBranch.GetCellData().AddArray(reorderedATPArray)

        print "There are", reorderedATPBranch.GetNumberOfPoints(), "ATP points for label", label, "..."
        print "There are", reorderedATPBranch.GetNumberOfCells(), "ATP cells for label", label, "..."

        # Write the VTK EC mesh file.
        reorderedMeshWriter = vtk.vtkXMLPolyDataWriter()
        reorderedMeshWriter.SetInput(reorderedATPBranch)
        reorderedMeshWriter.SetFileName(atpVTKFiles[label])
        reorderedMeshWriter.Update()

    parentFile.close()
    leftBranchFile.close()
    rightBranchFile.close()

    print "All done ..."

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)