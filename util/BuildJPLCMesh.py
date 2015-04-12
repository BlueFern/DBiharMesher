# -*- coding: utf-8 -*-
"""
Create an initial JPLC Profile for an ECs Mesh and write it out as .vtp.
"""

# TODO: Separating the dataset into branches would make the script simpler.
# Branch labelling is also required for processing tree datasets of arbitrary complexity.

import os
import vtk
import numpy
import matplotlib.pyplot as pyplot

# This is for the c216 mesh.
'''
os.chdir("/home/cza14/BlueFern/WrapDbihar/tmpData/c216")
meshFile = "quadMeshFullECc216.vtp"
jplcFile = "quadMeshFullJPLCc216.vtp"
numBranches = 3
numQuads = 216
numECsPerCol = 4
jplcGradient = 0.15
''' and None

# This is for the c4032 mesh.
'''
os.chdir("/home/cza14/BlueFern/WrapDbihar/tmpData/c4032")
meshFile = "quadMeshFullECc4032.vtp"
jplcFile = "quadMeshFullJPLCc4032.vtp"
numBranches = 3
numQuads = 4032
numECsPerCol = 4
jplcGradient = 0.05
''' and None

# This is for the c4080 mesh.
# '''
os.chdir("/home/cza14/BlueFern/WrapDbihar/tmpData/c4080")
meshFile = "quadMeshFullECc4080.vtp"
jplcFile = "quadMeshFullJPLCc4080.vtp"
numBranches = 3
numQuads = 4080
numECsPerCol = 4
jplcGradient = 0.03
# ''' and None

jplcMin = 0.2
jplcMax = 2.7

# Sigmoid function for providing JPLC values. The jplcGradient variable
# controls the "spread" of the values across the given domain.
def sigmoidJPLC(x):
    return jplcMin + (jplcMax / (1.0 + numpy.exp(-jplcGradient * x)))

def main():
    # Report our CWD just for testing purposes.
    print "CWD:", os.getcwd()    

    # Read in the mesh.
    meshReader = vtk.vtkXMLPolyDataReader()
    meshReader.SetFileName(meshFile)
    meshReader.Update()

    mesh = meshReader.GetOutput()
    
    # Put it through centroids filter.
    # Use VTK centroid filter to get the centroids in the right order
    # from the reorderedECMeshBranch.
    centroidFilter = vtk.vtkCellCenters()
    centroidFilter.SetInput(mesh)

    # Create a vertex cell for each point.
    pointsToVerticesFilter = vtk.vtkVertexGlyphFilter()
    pointsToVerticesFilter.SetInputConnection(centroidFilter.GetOutputPort())
    pointsToVerticesFilter.Update()
    
    jplcDataset = pointsToVerticesFilter.GetOutput()
    
    # Here we are assuming that cell ordering has been preserved.
    gridCoords = mesh.GetCellData().GetArray("gridCoords")

    # Create a new cell data array "parametricDistance" in the centroids
    # dataset with the values from the previous step. This array is really not
    # needed for any computations is good to have for verification. Also,
    # remember the max value for each branch.
    axialDist = vtk.vtkIntArray()
    axialDist.SetName("parametriDistance")
    
    # Iterate over all quad cells in the mesh and pull out the first component
    # of the "gridCoords" cell data array. Here we assuming the mesh contains
    # only quad cells.
    # If the branches were labelled properly here, it would be much simpler
    # to set the neganive axiial distance values by examining the range first.
    for cellId in range(0, mesh.GetNumberOfCells()):
        axialDist.InsertNextValue(gridCoords.GetTuple(cellId)[0])

    jplcArray = vtk.vtkDoubleArray()
    jplcArray.SetName("initialJPLC")
    
    axialDistRange = [0, 0]
    axialDist.GetRange(axialDistRange, 0)
    
    # Iterate over the "axialDist" cell array and for each value call
    # the agonistValue function and store the computed values in the new
    # "initialJPLC" cell data array.
    for cellId in range(0, axialDist.GetNumberOfTuples()):
        # Figure out the current branch id.
        branchId = cellId / (axialDist.GetNumberOfTuples() / 3)
        
        # For the first branch the parametric distance value must be in the
        # range [maxDist - 1, 0). Range values are required for this.
        
        # TODO: The JPLC calculations here are repeated many times for the same value
        # of the parametric distnace. It would make sense to have a lookup table for previously
        # calculated values.
        if branchId == 0:
            distVal = axialDist.GetValue(cellId) - axialDistRange[1] - 1
            jplcVal = sigmoidJPLC(distVal)
            jplcArray.InsertNextValue(jplcVal)
        else:
            distVal = axialDist.GetValue(cellId)
            jplcVal = sigmoidJPLC(distVal)
            jplcArray.InsertNextValue(jplcVal)

    # Assert the number of cells is equal to the number of items in the cell arrays.
    assert axialDist.GetNumberOfTuples() == jplcDataset.GetNumberOfCells(), "Number of cells (%d) and cell data values (%d) mismatch." % (axialDist.GetNumberOfTuples(), jplcDataset.GetNumberOfCells)
    assert jplcArray.GetNumberOfTuples() == jplcDataset.GetNumberOfCells(), "Number of cells (%d) and cell data values (%d) mismatch." % (jplcArray.GetNumberOfTuples(), jplcDataset.GetNumberOfCells())
    
    jplcDataset.GetCellData().AddArray(axialDist)
    jplcDataset.GetCellData().AddArray(jplcArray)
    
    jplcMapWriter = vtk.vtkXMLPolyDataWriter()
    jplcMapWriter.SetFileName(jplcFile)
    jplcMapWriter.SetInput(jplcDataset)
    jplcMapWriter.Update()
    
    # Provide a quick visualisation of the JPLC profile for validation.
    pointsX = range(-int(axialDistRange[1] - 1), int(axialDistRange[1]))
    pointsY = []
    for pt in pointsX:
        pointsY.append(sigmoidJPLC(pt))
    
    pyplot.plot(pointsX, pointsY, 'b')
    pyplot.show()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)