# -*- coding: utf-8 -*-
"""
Create an initial V2 ATP Profile for an ECs Mesh and write it out as .vtp.
"""
import os
import sys
import math
import vtk
import numpy
import matplotlib.pyplot as pyplot

centrelineFile = None
meshFile = None
debugAtpFile = None
atpFile = None
numBranches = 3
numQuads = 0
numAxialQuads = 0
numECsPerCol = 0
atpGradient = 0
atpMin = 0.1
atpMax = 1.0
outMin = -1.0
outMax = 1.0

def rescale(val, inMin, inMax):
    return (val - inMin) * (outMax - outMin) / (inMax - inMin) + outMin

# Sigmoid function for providing ATP values. The atpGradient variable
# controls the "spread" of the values across the given domain.
def sigmoidATP(x):
    return atpMin + (atpMax / (1.0 + numpy.exp(-atpGradient * x)))

def buildATPMesh():
    # Report our CWD just for testing purposes.
    print "CWD:", os.getcwd()

    # Read in the mesh.
    print 'Reading', meshFile
    meshReader = vtk.vtkXMLPolyDataReader()
    meshReader.SetFileName(meshFile)
    meshReader.Update()

    ecMesh = meshReader.GetOutput()

    # Read in the centreline.
    print 'Reading', centrelineFile
    centrelineReader = vtk.vtkPolyDataReader()
    centrelineReader.SetFileName(centrelineFile)
    centrelineReader.Update()

    centreline = centrelineReader.GetOutput()
    origin = centreline.GetPoint(0)
    print 'origin:', origin

    # Put the ecMesh through centroids filter.
    centroidFilter = vtk.vtkCellCenters()
    centroidFilter.SetInput(ecMesh)
    centroidFilter.Update()

    centroids = centroidFilter.GetOutput()

    # Iterate over each centroid and find the closest segment
    centroidPoints = centroids.GetPoints()

    # Only for DEBUG output.
    distArray = vtk.vtkDoubleArray()
    distArray.SetName("Dist")

    # For each point calculate the distance from origin.
    totalPoints = centroidPoints.GetNumberOfPoints()
    for ptId in range(totalPoints):        
        distance = vtk.vtkMath.Distance2BetweenPoints(origin, centroidPoints.GetPoint(ptId))
        distArray.InsertNextValue(math.sqrt(distance))

    # Get the range of the distance values.
    inMin, inMax = distArray.GetRange()
    
    atpArray = vtk.vtkFloatArray()
    atpArray.SetName('initialATP')
    
    # Normalise distance values.
    for i in range(distArray.GetNumberOfTuples()):
        dist = distArray.GetTuple(i)[0]
        distRescaled = rescale(dist, inMin, inMax)
        atpVal = sigmoidATP(distRescaled)
        atpArray.InsertNextValue(atpVal)
    
    # Prepare debug ATP mesh.
    debugAtpDataset = ecMesh
    debugAtpDataset.GetCellData().AddArray(distArray)
    debugAtpDataset.GetCellData().AddArray(atpArray)

    # Save the debug ATP mesh.
    debugAtpMapWriter = vtk.vtkXMLPolyDataWriter()
    debugAtpMapWriter.SetFileName(debugAtpFile)
    debugAtpMapWriter.SetInput(debugAtpDataset)
    debugAtpMapWriter.Update()

    # Prepare the ATP mesh by converting all points to vercices.
    pointsToVerticesFilter = vtk.vtkVertexGlyphFilter()
    pointsToVerticesFilter.SetInput(centroids)
    pointsToVerticesFilter.Update()
    atpDataset = pointsToVerticesFilter.GetOutput()
    atpDataset.GetCellData().AddArray(atpArray)

    # Assert the number of cells is equal to the number of items in the cell arrays.
    assert atpArray.GetNumberOfTuples() == debugAtpDataset.GetNumberOfCells(), "Number of cells (%d) and cell data values (%d) mismatch." % (atpArray.GetNumberOfTuples(), debugAtpDataset.GetNumberOfCells())

    atpMapWriter = vtk.vtkXMLPolyDataWriter()
    atpMapWriter.SetFileName(atpFile)
    atpMapWriter.SetInput(atpDataset)
    atpMapWriter.Update()
    
    # Provide a quick visualisation of the ATP profile for validation.
    pointsX = numpy.arange(outMin, outMax, 0.001)
    pointsY = []
    for pt in pointsX:
        pointsY.append(sigmoidATP(pt))
    
    pyplot.plot(pointsX, pointsY, 'b')
    pyplot.show()

def usage():
    print 'This script is to be run with global parameters (input centrelin, EC mesh, etc.) set in the calling script.'

if __name__ == '__main__':
    print 'Starting', os.path.basename(__file__)
    usage()
    print 'Exiting', os.path.basename(__file__)