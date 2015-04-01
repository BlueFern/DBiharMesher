# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 17:09:57 2015
"""
import os

import vtk

meshSet0 = [
"quadMeshFullc216.vtp",
"quadMeshFullECc216.vtp",
"quadMeshFullSMCc216.vtp"
]
numQuadsPerRing0 = 12
numECsPerQuad = 4
numSMCsPerQuad = 4

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

taskTXTFiles = [
"txt/parent_points.txt",
"txt/parent_cells.txt",
"txt/left_daughter_points.txt",
"txt/left_daughter_cells.txt",
"txt/right_daughter_points.txt",
"txt/right_daughter_cells.txt"
]

ecCentroidTXTFiles = [
"txt/parent_ec_centeroid_points.txt",
"txt/parent_ec_centeroid_cells.txt",
"txt/left_daughter_ec_centeroid_points.txt",
"txt/left_daughter_ec_centeroid_cells.txt",
"txt/right_daughter_ec_centeroid_points.txt",
"txt/right_daughter_ec_centeroid_cells.txt"
]

ecTXTFiles = [
"txt/parent_ec_mesh_points.txt",
"txt/parent_ec_mesh_cells.txt",
"txt/left_daughter_ec_mesh_points.txt",
"txt/left_daughter_ec_mesh_cells.txt",
"txt/right_daughter_ec_mesh_points.txt",
"txt/right_daughter_ec_mesh_cells.txt"
]

smcTXTFiles = [
"txt/parent_smc_mesh_points.txt",
"txt/parent_smc_mesh_cells.txt",
"txt/left_daughter_smc_mesh_points.txt",
"txt/left_daughter_smc_mesh_cells.txt",
"txt/right_daughter_smc_mesh_points.txt",
"txt/right_daughter_smc_mesh_cells.txt"
]

def main():
    # This is where the data is for testing purposes.
    os.chdir("/home/cza14/BlueFern/LocalData/c216")
    print "Current working directory:", os.getcwd()    

    # Working with the task mesh.
    # Working with the task mesh.
    # Working with the task mesh.
    taskMeshReader = vtk.vtkXMLPolyDataReader()
    taskMeshReader.SetFileName(meshSet0[0])
    taskMeshReader.Update()

    taskMesh = taskMeshReader.GetOutput()
    
    # Get the range of branch labels.
    labelRange = [0, 0]
    taskMesh.GetCellData().GetScalars().GetRange(labelRange, 0)
    
    # Convert label range to a list of labels.
    labelRange = range(int(labelRange[0]), int(labelRange[1]) + 1)    
    print "Labels found in task mesh:", labelRange
   
    parentPointsFile = open(taskTXTFiles[0], 'w')
    parentCellsFile = open(taskTXTFiles[1], 'w')

    leftPointsFile = open(taskTXTFiles[2], 'w')
    leftCellsFile = open(taskTXTFiles[3], 'w')
    
    rightPointsFile = open(taskTXTFiles[4], 'w')
    rightCellsFile = open(taskTXTFiles[5], 'w')
   
    # For every label in the range of labels we want to extract all cells.
    for label in labelRange:
        # Use this filter to extract the cells for a given label value.
        branchSelector = vtk.vtkThreshold()
        branchSelector.SetInput(taskMesh)
        branchSelector.ThresholdBetween(label,label);
        branchSelector.Update()

        taskMeshBranch = branchSelector.GetOutput()
        
        print "There are", taskMeshBranch.GetNumberOfCells(), \
        "cells for label", label, "..."
        
        # Create new vtkPolyData object for the new reordered mesh.
        reorderedTaskMeshBranch = vtk.vtkPolyData()
        # Use the old points.
        reorderedTaskMeshBranch.SetPoints(taskMeshBranch.GetPoints())
        
        # New vtkCellArray for storing reordeced cells.
        reorderedCellArray = vtk.vtkCellArray()
        numRings = taskMeshBranch.GetNumberOfCells() / numQuadsPerRing0;
        ringIds = range(0, numRings);
        ringIds.reverse()

        # Decide which files to write to.
        pointsOf = ''
        cellsOf = ''
        if label == 0:
            pointsOf = parentPointsFile
            cellsOf = parentCellsFile
        elif label == 1:
            pointsOf = leftPointsFile
            cellsOf = leftCellsFile
        elif label == 2:
            pointsOf = rightPointsFile
            cellsOf = rightCellsFile
       
        # Iterante over the rings in reverse order.
        for ringNum in ringIds:
            # Iterate over the cells in normal order.
            for cellNum in range(0, numQuadsPerRing0):
                # Calculate the 'real' cell id and the corresponding cell.
                cellId = ringNum * numQuadsPerRing0 + cellNum
                cell = taskMeshBranch.GetCell(cellId)
                reorderedCellArray.InsertNextCell(cell)
                
                pointIdList = [cell.GetNumberOfPoints()]
                for pPos in range(0, cell.GetNumberOfPoints()):
                    pointIdList.append(cell.GetPointId(pPos))
                    
                pointIdListStr = ' '.join(str(i) for i in pointIdList)
                cellsOf.write(pointIdListStr + '\n')
                
                # Write points to text file.
                for pPos in range(0, cell.GetNumberOfPoints()):
                    writePoint = False
                    if ringNum == ringIds[0]:
                        if cellNum == 0:
                            writePoint = True
                        elif pPos == 1 or pPos == 2:
                            writePoint = True
                    else:
                        if cellNum == 0:
                            if pPos == 0 or pPos == 1:
                                writePoint = True
                        else:
                            if pPos == 1:
                                writePoint = True
                    if writePoint == True:
                        print pPos,
                        point = taskMeshBranch.GetPoint(cell.GetPointId(pPos))
                        pointStr = ' '.join(format(i , '.6f') for i in point)
                        pointsOf.write(pointStr + '\n')
            print '\n'

        # Put the reordered cells into the reordered mesh.
        reorderedTaskMeshBranch.SetPolys(reorderedCellArray)
        
        # Write the VTK file.
        reorderedMeshWriter = vtk.vtkXMLPolyDataWriter()
        reorderedMeshWriter.SetInput(reorderedTaskMeshBranch)
        reorderedMeshWriter.SetFileName(taskVTKFiles[label])
        reorderedMeshWriter.Update()

    parentPointsFile.close()
    parentCellsFile.close()
    
    leftPointsFile.close()
    leftCellsFile.close()
    
    rightPointsFile.close()
    rightCellsFile.close()

    # Working with EC mesh.
    # Working with EC mesh.
    # Working with EC mesh.
    ecMeshReader = vtk.vtkXMLPolyDataReader()
    ecMeshReader.SetFileName(meshSet0[1])
    ecMeshReader.Update()

    ecMesh = ecMeshReader.GetOutput()
 
    
    # Working with SMC mesh.
    # Working with SMC mesh.
    # Working with SMC mesh.
    smcMeshReader = vtk.vtkXMLPolyDataReader()
    smcMeshReader.SetFileName(meshSet0[2])
    smcMeshReader.Update()

    smcMesh = smcMeshReader.GetOutput()
    
if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)