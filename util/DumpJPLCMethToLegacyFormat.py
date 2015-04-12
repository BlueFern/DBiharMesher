# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 17:09:57 2015

This script read an ECs mesh produced by WrapDbihar code and converts it to
a JPLC profile map for the legacy Coupled Cells code.

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
os.chdir('/home/cza14/BlueFern/WrapDbihar/tmpData/c216')
numQuadsPerRing0 = 12
taskMeshIn = "quadMeshFullc216.vtp"
jplcMeshIn = "quadMeshFullJPLCc216.vtp"
''' and None

'''
#os.chdir('/home/cza14/BlueFern/WrapDbihar/tmpData/c4032')
numQuadsPerRing0 = 48
taskMeshIn = "quadMeshFullc4032.vtp"
jplcMeshIn = "quadMeshFullJPLCc4032.vtp"
''' and None

#'''
os.chdir('/home/cza14/BlueFern/WrapDbihar/tmpData/c4080')
numQuadsPerRing0 = 40
taskMeshIn = "quadMeshFullc4080.vtp"
jplcMeshIn = "quadMeshFullJPLCc4080.vtp"
#''' and None

'''
#os.chdir('/home/cza14/BlueFern/WrapDbihar/tmpData/c8112')
numQuadsPerRing0 = 40
taskMeshIn = "quadMeshFullc4080.vtp"
jplcMeshIn = "quadMeshFullJPLCc8112.vtp"
''' and None


# VTK files to write.
jplcVTKFiles = [
"vtk/jplc_parent.vtp",
"vtk/jplc_left_daughter.vtp",
"vtk/jplc_right_daughter.vtp",
]

jplcTXTFiles = [
"txt/parent_jplc.txt",
"txt/left_daughter_jplc.txt",
"txt/right_daughter_jplc.txt",
]

def main():
    # This is where the data is for testing purposes.
    print "Current working directory:", os.getcwd()

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
    ecMeshReader = vtk.vtkXMLPolyDataReader()
    ecMeshReader.SetFileName(jplcMeshIn)
    ecMeshReader.Update()

    # Original ECs mesh to work with.
    jplcMesh = ecMeshReader.GetOutput()
    print "There are", jplcMesh.GetNumberOfCells(), "JPLC values in total ..."

    # Prepare the stupid TXT files for output.    
    parentFile = open(jplcTXTFiles[0], 'w')
    leftBranchFile = open(jplcTXTFiles[1], 'w')
    rightBranchFile = open(jplcTXTFiles[2], 'w')

    # For every label in the range of labels we want to extract all ECs.
    for label in labelRange:

        # Keep track of how many branches we need to skip.
        numECsPerLabel = numQuadsPerRing0 * numRingsPerLabel[label] * numECsPerQuad
        jplcCellOffset = label * numECsPerLabel

        print "jplcCellOffset", jplcCellOffset

        # Collect cell ids to select.
        selectionIds = vtk.vtkIdTypeArray()
        for sId in range(0, numECsPerLabel):
            selectionIds.InsertNextValue(jplcCellOffset + sId)

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
        selectionExtractor.SetInput(0, jplcMesh)
        selectionExtractor.SetInput(1, selection)
        selectionExtractor.Update()

        extractedCells = selectionExtractor.GetOutput()

        # Ring ids list for traversal.
        ringIds = range(0, numRingsPerLabel[label])
        ringIds.reverse()

        # Number of ECs rows is the number of ECs per quad.
        rowIds = range(0, numECsPerCol)
        rowIds.reverse()

        # New vtkCellArray for storing reordeced cells.
        reorderedCellArray = vtk.vtkCellArray()
        
        reorderedJPLCArray = vtk.vtkDoubleArray()
        reorderedJPLCArray.SetName("initialJPLC")

        # Decide which TXT files to write to.
        pointsOf = ''
        
        if label == 0:
            pointsOf = parentFile
        elif label == 1:
            pointsOf = leftBranchFile
        elif label == 2:
            pointsOf = rightBranchFile

        print "Writing TXT file for ECs JPLC:"
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
                        cell = extractedCells.GetCell(realId)
                        reorderedCellArray.InsertNextCell(cell)
                        
                        jplcVal = extractedCells.GetCellData().GetArray("initialJPLC").GetValue(realId)
                        
                        reorderedJPLCArray.InsertNextValue(jplcVal)

                        # Write the value to the txt file.
                        pointsOf.write(format(jplcVal, '.6f') + '\n')
                        
        # Create new vtkPolyData object for the new reordered mesh.
        reorderedJPLCBranch = vtk.vtkPolyData()

        # Insert our new points.
        reorderedJPLCBranch.SetPoints(extractedCells.GetPoints())

        # Set the reordered cells to the reordered JPLC mesh.
        reorderedJPLCBranch.SetVerts(reorderedCellArray)
        
        # Set the reordered JPLC values to the reordered JPLC mesh.
        reorderedJPLCBranch.GetCellData().AddArray(reorderedJPLCArray)

        print "There are", reorderedJPLCBranch.GetNumberOfPoints(), "JPLC points for label", label, "..."
        print "There are", reorderedJPLCBranch.GetNumberOfCells(), "JPLC cells for label", label, "..."

        # Write the VTK EC mesh file.
        reorderedMeshWriter = vtk.vtkXMLPolyDataWriter()
        reorderedMeshWriter.SetInput(reorderedJPLCBranch)
        reorderedMeshWriter.SetFileName(jplcVTKFiles[label])
        reorderedMeshWriter.Update()

    parentFile.close()
    leftBranchFile.close()
    rightBranchFile.close()

    print "All done ..."

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)