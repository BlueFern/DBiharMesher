# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 17:09:57 2015

This script converts output from WrapDbihar code to the format expected by the
legacy Coupled Cells code.

The code reorders cells and produces the required files in TXT format.
VTK files are written out for visual verification.

This code assumes there are three branches (parent and two daughters) in the data.

"""

import os
import vtk

numECsPerCol = 4
numSMCsPerRow = 4

numECsPerRow = numSMCsPerRow * 5
numSMCsPerCol = numECsPerCol * 13

numECsPerQuad = numECsPerRow * numECsPerCol
numSMCsPerQuad = numSMCsPerCol * numSMCsPerRow

numQuadsPerRing = 0
meshSet = []

# VTK files to write.
taskVTKFiles = [
"vtk/parent.vtp",
"vtk/left_daughter.vtp",
"vtk/right_daughter.vtp",
]

ecCentroidVTKFiles = [
"vtk/ec_centeroid_parent.vtp",
"vtk/ec_centeroid_left_daughter.vtp",
"vtk/ec_centeroid_right_daughter.vtp",
]

ecVTKFiles = [
"vtk/ec_mesh_parent.vtp",
"vtk/ec_mesh_left_daughter.vtp",
"vtk/ec_mesh_right_daughter.vtp",
]

smcVTKFiles = [
"vtk/smc_mesh_parent.vtp",
"vtk/smc_mesh_left_daughter.vtp",
"vtk/smc_mesh_right_daughter.vtp"
]

def writeLegacyVTK():
    # This is where the data is for testing purposes.
    print "Current working directory:", os.getcwd()
    
    if os.path.isdir("vtk") == False:
        os.makedirs("vtk")
        print "Cretated vtk output directory..."
    
    if os.path.isdir("files") == False:
        os.makedirs("files")
        print "Created files ouptut directory..."
        

    # Working with the task mesh.
    taskMeshReader = vtk.vtkXMLPolyDataReader()
    taskMeshReader.SetFileName(meshSet[0])
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
        branchSelector.SetInputData(taskMesh)
        branchSelector.ThresholdBetween(label,label);
        branchSelector.Update()

        taskMeshBranch = branchSelector.GetOutput()

        # New vtkPoints for storing reordered points.
        reorderedPoints = vtk.vtkPoints()

        # New vtkCellArray for storing reordeced cells.
        reorderedCellArray = vtk.vtkCellArray()
        numQuadRowsPerBranch = taskMeshBranch.GetNumberOfCells() / numQuadsPerRing;
        numRingsPerLabel[label] = numQuadRowsPerBranch
        ringIds = range(0, numQuadRowsPerBranch);

        # Working with rows in reverse order: UPSTREAM.
        ringIds.reverse()


        rowBase = 0
        # Iterate over the rings in reverse order.
        for ringNum in ringIds:
            # Iterate over the cells in normal order.
            for cellNum in range(0, numQuadsPerRing):
                # Calculate the 'real' cell id and get the corresponding cell.
                cellId = ringNum * numQuadsPerRing + cellNum
                cell = taskMeshBranch.GetCell(cellId)

                # The ids to be written to the TXT file.
                pointIdList = [cell.GetNumberOfPoints()]

                # Write the appropriate points to TXT file.
                for pPos in range(0, cell.GetNumberOfPoints()):
                    newPoint = False
                    if ringNum == ringIds[0]:
                        if cellNum == 0:
                            newPoint = True
                        elif pPos == 1 or pPos == 2:
                            newPoint = True
                    else:
                        if cellNum == 0:
                            if pPos == 0 or pPos == 1:
                                newPoint = True
                        else:
                            if pPos == 1:
                                newPoint = True

                    if newPoint == True:

                        # Inserting a new point...
                        point = taskMeshBranch.GetPoint(cell.GetPointId(pPos))
                        # ... with a new id.
                        newId = reorderedPoints.InsertNextPoint(point)
                        pointIdList.append(newId)

                        # To make it easier for remembering the number of points instered in a row.
                        if cellNum == 0 and pPos == 0:
                            rowBasePrev = newId
                    else:
                        # Perhaps this can be done in a nicer way.
                        # Calculate the id of a previously inserted point.
                        if ringNum == ringIds[0]:
                            if cellNum == 1:
                                if pPos == 0:
                                    pointIdList.append(1L)
                                elif pPos == 3:
                                    pointIdList.append(2L)
                            else:
                                if pPos == 0:
                                    pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 2))
                                elif pPos == 3:
                                    pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 3))
                        elif ringNum == ringIds[1]:
                            if cellNum == 0:
                                if pPos == 2:
                                    pointIdList.append(1L)
                                elif pPos == 3:
                                    pointIdList.append(0L)
                            else:
                                if pPos == 0:
                                    pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 1))
                                elif pPos == 2:
                                    pointIdList.append(long(cellNum * 2 + 2))
                                elif pPos == 3:
                                    if cellNum == 1:
                                        pointIdList.append(1L)
                                    else:
                                        pointIdList.append(long(cellNum * 2))
                        else:
                            if cellNum == 0:
                                if pPos == 2:
                                    pointIdList.append(long(rowBase + 1))
                                elif pPos == 3:
                                    pointIdList.append(long(rowBase))
                            else:
                                if pPos == 0:
                                    pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 1))
                                elif pPos == 2:
                                    pointIdList.append(long(rowBase + cellNum + 1))
                                elif pPos == 3:
                                    pointIdList.append(long(rowBase + cellNum))

                # print pointIdList, rowBase

                # Insert the ids into the cell array.
                newCell = vtk.vtkQuad()
                newCell.GetPointIds().Reset()
                for id in pointIdList[1:]:
                    newCell.GetPointIds().InsertNextId(id)
                reorderedCellArray.InsertNextCell(newCell)

            rowBase = rowBasePrev

        # print '\n'
        print "Inserted", reorderedPoints.GetNumberOfPoints(), "task mesh points for label", label, "..."
        print "Inserted", reorderedCellArray.GetNumberOfCells(), "task mesh cells for label", label, "..."

        # Create new vtkPolyData object for the new reordered mesh.
        reorderedTaskMeshBranch = vtk.vtkPolyData()

        # Put the reordered points and cells into the reordered mesh.
        reorderedTaskMeshBranch.SetPoints(reorderedPoints)
        reorderedTaskMeshBranch.SetPolys(reorderedCellArray)

        # Write the VTK file.
        reorderedMeshWriter = vtk.vtkXMLPolyDataWriter()
        reorderedMeshWriter.SetInputData(reorderedTaskMeshBranch)
        reorderedMeshWriter.SetFileName(taskVTKFiles[label])
        reorderedMeshWriter.Update()

    print "Rings per label:", numRingsPerLabel, "..."
    ringsPerLabelVals = numRingsPerLabel.values()

    # Check all rings per label values are the same.
    assert ringsPerLabelVals[1:] == ringsPerLabelVals[:-1], "All values of rings per label must be identical. Generated output is invalid ..."

    # Working with EC mesh.

    ecMeshReader = vtk.vtkXMLPolyDataReader()
    ecMeshReader.SetFileName(meshSet[1])
    ecMeshReader.Update()

    # Original ECs mesh to work with.
    ecMesh = ecMeshReader.GetOutput()
    print "There are", ecMesh.GetNumberOfCells(), "ECs in total ..."

    # For every label in the range of labels we want to extract all ECs.
    for label in labelRange:

        # Keep track of how many branches we need to skip.
        numECsPerLabel = numQuadsPerRing * numRingsPerLabel[label] * numECsPerQuad
        ecCellOffset = label * numECsPerLabel

        print "ecCellOffset", ecCellOffset

        # Collect cell ids to select.
        selectionIds = vtk.vtkIdTypeArray()
        for sId in range(0, numECsPerLabel):
            selectionIds.InsertNextValue(ecCellOffset + sId)

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
        selectionExtractor.SetInputData(0, ecMesh)
        selectionExtractor.SetInputData(1, selection)
        selectionExtractor.Update()

        extractedECs = selectionExtractor.GetOutput()

        # Ring ids list for traversal.
        ringIds = range(0, numRingsPerLabel[label])
        ringIds.reverse()

        # Number of ECs rows is the number of ECs per quad.
        rowIds = range(0, numECsPerCol)
        rowIds.reverse()

        # The ECs are organised in rings of blocks of cells.
        # New vtkCellArray for storing reordeced cells.
        reorderedCellArray = vtk.vtkCellArray()

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
                    for ecNum in range(0, numECsPerRow):
                        # Calculate the 'real' ec cell id and get the corresponding cell.
                        ecId = quadId * numECsPerQuad + rowNum * numECsPerRow + ecNum
                        ecCell = extractedECs.GetCell(ecId)
                        reorderedCellArray.InsertNextCell(ecCell)

        # Create new vtkPolyData object for the new reordered mesh.
        reorderedECMeshBranch = vtk.vtkPolyData()

        # Insert our new points.
        reorderedECMeshBranch.SetPoints(extractedECs.GetPoints())

        # Set the reordered cells to the reordered ECs mesh.
        reorderedECMeshBranch.SetPolys(reorderedCellArray)

        # New vtkPoints for storing reordered points.
        reorderedPoints = vtk.vtkPoints()

        # New vtkCellArray for storing reordeced cells.
        reorderedCellArray = vtk.vtkCellArray()

        rowBase = 0
        # Iterate over quads in normal order because they have been reordered.
        for quadNum in range(0, numRingsPerLabel[label] * numQuadsPerRing):
            # Iterate over rows in normal order because they have been reordered.
            for rowNum in range(0, numECsPerCol):
                # Iterate over the ECs in the row in normal order.
                for ecNum in range(0, numECsPerRow):
                    # Calculate the 'real' ec cell id and get the corresponding cell.
                    ecId = quadNum * numECsPerQuad + rowNum * numECsPerRow + ecNum
                    ecCell = reorderedECMeshBranch.GetCell(ecId)

                    # The ids to be written to the TXT file.
                    pointIdList = [ecCell.GetNumberOfPoints()]

                    # Write the appropriate points to the TXT file.
                    for pPos in range(0, ecCell.GetNumberOfPoints()):
                        newPoint = False
                        if rowNum == 0:
                            if ecNum == 0:
                                newPoint = True
                            elif pPos == 1 or pPos == 2:
                                newPoint = True
                        else:
                            if ecNum == 0:
                                if pPos == 0 or pPos == 1:
                                    newPoint = True
                            else:
                                if pPos == 1:
                                    newPoint = True

                        if newPoint == True:

                            # Inserting a new point...
                            point = reorderedECMeshBranch.GetPoint(ecCell.GetPointId(pPos))
                            # ... with a new id.
                            newId = reorderedPoints.InsertNextPoint(point)
                            pointIdList.append(newId)


                            if ecNum == 0 and pPos == 0:
                                rowBasePrev = newId
                        else:
                            # Perhaps this can be done in a nicer way.
                            # Calculate the ide of a previously inserted point.
                            if rowNum == 0:
                                if ecNum == 1:
                                    if pPos == 0:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 3))
                                    elif pPos == 3:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 4))
                                else:
                                    if pPos == 0:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 2))
                                    elif pPos == 3:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 3))
                            elif rowNum == 1:
                                if ecNum == 0:
                                    if pPos == 2:
                                        pointIdList.append(long(rowBase + 1))
                                    elif pPos == 3:
                                        pointIdList.append(long(rowBase))
                                else:
                                    if pPos == 0:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 1))
                                    elif pPos == 2:
                                        pointIdList.append(long(rowBase + ecNum * 2 + 2))
                                    elif pPos == 3:
                                        if ecNum == 1:
                                            pointIdList.append(long(rowBase + 1))
                                        else:
                                            pointIdList.append(long(rowBase + ecNum * 2))
                            else:
                                if ecNum == 0:
                                    if pPos == 2:
                                        pointIdList.append(long(rowBase + 1))
                                    elif pPos == 3:
                                        pointIdList.append(long(rowBase))
                                else:
                                    if pPos == 0:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 1))
                                    elif pPos == 2:
                                        pointIdList.append(long(rowBase + ecNum + 1))
                                    elif pPos == 3:
                                        pointIdList.append(long(rowBase + ecNum))

                    # print pointIdList, rowBase

                    # Insert the ids into the cell array.
                    newCell = vtk.vtkQuad()
                    newCell.GetPointIds().Reset()
                    for id in pointIdList[1:]:
                        newCell.GetPointIds().InsertNextId(id)
                    reorderedCellArray.InsertNextCell(newCell)

                rowBase = rowBasePrev

        # print '\n'
        print "There are", reorderedPoints.GetNumberOfPoints(), "ECs points for label", label, "..."
        print "There are", reorderedCellArray.GetNumberOfCells(), "ECs cells for label", label, "..."

        # Create new vtkPolyData object for the new reordered mesh.
        reorderedECs = vtk.vtkPolyData()

        # Put the reordered points and cells into the mesh.
        reorderedECs.SetPoints(reorderedPoints)
        reorderedECs.SetPolys(reorderedCellArray)

        # Write the VTK EC mesh file.
        reorderedMeshWriter = vtk.vtkXMLPolyDataWriter()
        reorderedMeshWriter.SetInputData(reorderedECs)
        reorderedMeshWriter.SetFileName(ecVTKFiles[label])
        reorderedMeshWriter.Update()

        # Use VTK centroid filter to get the centroids in the right order
        # from the reorderedECMeshBranch.
        centroidFilter = vtk.vtkCellCenters()
        centroidFilter.SetInputData(reorderedECs)
        centroidFilter.Update()

        # Create a vertex cell for each point.
        pointsToVerticesFilter = vtk.vtkVertexGlyphFilter()
        pointsToVerticesFilter.SetInputData(centroidFilter.GetOutput())
        pointsToVerticesFilter.Update()

        reorderedCentroidBranch = pointsToVerticesFilter.GetOutput()

        # Write the VTK EC centrouid file.
        centroidWriter = vtk.vtkXMLPolyDataWriter()
        centroidWriter.SetInputData(reorderedCentroidBranch)
        centroidWriter.SetFileName(ecCentroidVTKFiles[label])
        centroidWriter.Update()

        # Write the centroids to the TXT points and cells files.
        for cId in range(0, reorderedCentroidBranch.GetNumberOfCells()):
            centCell = reorderedCentroidBranch.GetCell(cId)
            centIds = [centCell.GetNumberOfPoints()]

            # Write centroid ids.
            ptId = centCell.GetPointId(0)
            centIds.append(ptId)

            # Write centroid points.
            point = reorderedCentroidBranch.GetPoint(ptId)




    # Working with SMC mesh.
    # Working with SMC mesh.
    # Working with SMC mesh.
    smcMeshReader = vtk.vtkXMLPolyDataReader()
    smcMeshReader.SetFileName(meshSet[2])
    smcMeshReader.Update()

    # Original SMCs mesh to work with.
    smcMesh = smcMeshReader.GetOutput()
    print "There are", smcMesh.GetNumberOfCells(), "SMCs in total ..."

    # For every label in the range of labels we want to extract all SMCs.
    for label in labelRange:

        # Keep track of how many branches we need to skip.
        numSMCsPerLabel = numQuadsPerRing * numRingsPerLabel[label] * numSMCsPerQuad
        smcCellOffset = label * numSMCsPerLabel

        print "smcCellOffset", smcCellOffset

        # Collect cell ids to select.
        selectionIds = vtk.vtkIdTypeArray()
        for sId in range(0, numSMCsPerLabel):
            selectionIds.InsertNextValue(smcCellOffset + sId)

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
        selectionExtractor.SetInputData(0, smcMesh)
        selectionExtractor.SetInputData(1, selection)
        selectionExtractor.Update()

        extractedSMCs = selectionExtractor.GetOutput()

        # Ring ids list for traversal.
        ringIds = range(0, numRingsPerLabel[label])
        ringIds.reverse()

        # Number of SMCs rows is the number of ECs per quad times 13.
        rowIds = range(0, numSMCsPerCol)
        rowIds.reverse()

        # The SMCs are organised in rings of blocks of cells.
        # New vtkCellArray for storing reordeced cells.
        reorderedCellArray = vtk.vtkCellArray()

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
                    for smcNum in range(0, numSMCsPerRow):
                        # Calculate the 'real' smc cell id and get the corresponding cell.
                        smcId = quadId * numSMCsPerQuad + rowNum * numSMCsPerRow + smcNum
                        smcCell = extractedSMCs.GetCell(smcId)
                        reorderedCellArray.InsertNextCell(smcCell)

        # Create new vtkPolyData object for the new reordered mesh.
        reorderedSMCMeshBranch = vtk.vtkPolyData()

        # Insert our new points.
        reorderedSMCMeshBranch.SetPoints(extractedSMCs.GetPoints())

        # Set the reordered cells to the reordered SMCs mesh.
        reorderedSMCMeshBranch.SetPolys(reorderedCellArray)

        # New vtkPoints for storing reordered points.
        reorderedPoints = vtk.vtkPoints()

        # New vtkCellArray for storing reordeced cells.
        reorderedCellArray = vtk.vtkCellArray()

        rowBase = 0
        # Iterate over quads in normal order because they have been reordered.
        for quadNum in range(0, numRingsPerLabel[label] * numQuadsPerRing):
            # Iterate over rows in normal order because they have been reordered.
            for rowNum in range(0, numSMCsPerCol):
                # Iterate over the SMCs in the row in normal order.
                for smcNum in range(0, numSMCsPerRow):
                    # Calculate the 'real' smc cell id and get the corresponding cell.
                    smcId = quadNum * numSMCsPerQuad + rowNum * numSMCsPerRow + smcNum
                    smcCell = reorderedSMCMeshBranch.GetCell(smcId)

                    # The ids to be written to the TXT file.
                    pointIdList = [smcCell.GetNumberOfPoints()]

                    # Write the appropriate points to the TXT file.
                    for pPos in range(0, smcCell.GetNumberOfPoints()):
                        newPoint = False
                        if rowNum == 0:
                            if smcNum == 0:
                                newPoint = True
                            elif pPos == 1 or pPos == 2:
                                newPoint = True
                        else:
                            if smcNum == 0:
                                if pPos == 0 or pPos == 1:
                                    newPoint = True
                            else:
                                if pPos == 1:
                                    newPoint = True
                        if newPoint == True:

                            # Inserting a new point...
                            point = reorderedSMCMeshBranch.GetPoint(smcCell.GetPointId(pPos))
                            # with a new id.
                            newId = reorderedPoints.InsertNextPoint(point)
                            pointIdList.append(newId)

                            if smcNum == 0 and pPos == 0:
                                rowBasePrev = newId
                        else:
                            # Perhaps this can be done in a nicer way.
                            # Calculate the ide of a previously inserted point.
                            if rowNum == 0:
                                if smcNum == 1:
                                    if pPos == 0:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 3))
                                    elif pPos == 3:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 4))
                                else:
                                    if pPos == 0:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 2))
                                    elif pPos == 3:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 3))
                            elif rowNum == 1:
                                if smcNum == 0:
                                    if pPos == 2:
                                        pointIdList.append(long(rowBase + 1))
                                    elif pPos == 3:
                                        pointIdList.append(long(rowBase))
                                else:
                                    if pPos == 0:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 1))
                                    elif pPos == 2:
                                        pointIdList.append(long(rowBase + smcNum * 2 + 2))
                                    elif pPos == 3:
                                        if smcNum == 1:
                                            pointIdList.append(long(rowBase + 1))
                                        else:
                                            pointIdList.append(long(rowBase + smcNum * 2))
                            else:
                                if smcNum == 0:
                                    if pPos == 2:
                                        pointIdList.append(long(rowBase + 1))
                                    elif pPos == 3:
                                        pointIdList.append(long(rowBase))
                                else:
                                    if pPos == 0:
                                        pointIdList.append(long(reorderedPoints.GetNumberOfPoints() - 1))
                                    elif pPos == 2:
                                        pointIdList.append(long(rowBase + smcNum + 1))
                                    elif pPos == 3:
                                        pointIdList.append(long(rowBase + smcNum))

                    # print pointIdList, rowBase

                    # Insert the ids into the cell array.
                    newCell = vtk.vtkQuad()
                    newCell.GetPointIds().Reset()
                    for id in pointIdList[1:]:
                        newCell.GetPointIds().InsertNextId(id)
                    reorderedCellArray.InsertNextCell(newCell)

                rowBase = rowBasePrev

        # print '\n'
        print "There are", reorderedPoints.GetNumberOfPoints(), "SMCs points for label", label, "..."
        print "There are", reorderedCellArray.GetNumberOfCells(), "SMCs cells for label", label, "..."

        # Create new vtkPolyData object for the new reordered mesh.
        reorderedSMCs = vtk.vtkPolyData()

        # Put the reordered points and cells in to the mesh.
        reorderedSMCs.SetPoints(reorderedPoints)
        reorderedSMCs.SetPolys(reorderedCellArray)

        # Write the VTK SMC mesh file.
        reorderedMeshWriter = vtk.vtkXMLPolyDataWriter()
        reorderedMeshWriter.SetInputData(reorderedSMCs)
        reorderedMeshWriter.SetFileName(smcVTKFiles[label])
        reorderedMeshWriter.Update()

    print "All done ..."
    print "... Except the last configuration_info.txt file ..."

    configFile = open("files/configuration_info.txt", 'w')
    configFile.write("Processors information\n")
    configFile.write("Total number of points per branch (vtk points) = %d\t\tm = %d n = %d\n" \
    % ((numQuadsPerRing + 1) * (numRingsPerLabel[0] + 1), (numQuadsPerRing + 1), (numRingsPerLabel[0] + 1)))
    configFile.write("Total number of cells per branch (vtk cells) = %d\t\tm = %d n = %d\n" \
    % (numQuadsPerRing * numRingsPerLabel[0], numQuadsPerRing, numRingsPerLabel[0]))
    configFile.write("Total number of SMC mesh points per processor mesh (vtk points) = %d\t\tm = %d n = %d\n" \
    % ((numSMCsPerCol + 1) * (numSMCsPerRow + 1), (numSMCsPerCol + 1), (numSMCsPerRow + 1)))
    configFile.write("Total number of SMC mesh cells per processor mesh (vtk cells) = %d\t\tm = %d n = %d\n" \
    % (numSMCsPerCol * numSMCsPerRow, numSMCsPerCol, numSMCsPerRow))
    configFile.write("Total number of EC mesh points per processor mesh (vtk points) = %d\t\tm = %d n = %d\n" \
    % ((numECsPerCol + 1) * (numECsPerRow + 1), (numECsPerCol + 1), (numECsPerRow + 1)))
    configFile.write("Total number of EC mesh cells per processor mesh (vtk cells) = %d\t\tm = %d n = %d\n" \
    % (numECsPerCol *numECsPerRow, numECsPerCol, numECsPerRow))
    configFile.write("Total number of EC mesh centeroid points per processor mesh (vtk points) = %d\t\tm = %d n = %d\n" \
    % (numECsPerCol *numECsPerRow, numECsPerCol, numECsPerRow))
    configFile.write("Total number of EC mesh centeroid cells per processor mesh (vtk cells) = %d\t\tm = %d n = %d\n" \
    % (numECsPerCol *numECsPerRow, numECsPerCol, numECsPerRow))

    configFile.close()

    print "Now it is all done for real ..."

def main():
    print "This script is to be run with global parameters (input, output files, etc.) set in the calling script."

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)

