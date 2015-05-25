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

# Sigmoid function for providing ATP values. The atpGradient variable
# controls the "spread" of the values across the given domain.
def sigmoidATP(x):
    return atpMin + (atpMax / (1.0 + numpy.exp(-atpGradient * x)))

def getParametricDistance(resampledCentreline, point):

    # Iterate over all line segments in the resampled centreline and find the number of the closest segment.
    distance = sys.maxint
    pos = 0
    lineId = 0

    # Resample the centreline to EC resolution.
    lineArray = resampledCentreline.GetLines()
    lineArray.InitTraversal()
    line = vtk.vtkIdList()
    tmpLineId = 0

    lineLengths = {}
    t = None
    closestPoint = None

    while lineArray.GetNextCell(line):

        lineLengths[tmpLineId] = line.GetNumberOfIds() - 1

        # Iterate over each pair of adjacent points and find the distance to from the point to the current.
        for ptId in range(line.GetNumberOfIds() - 1):

            t = vtk.mutable(0)
            closestPoint = [0,0,0]
            # Current distance.
            tmpDistance = vtk.vtkLine.DistanceToLine(point, \
                                       resampledCentreline.GetPoint(line.GetId(ptId)), \
                                       resampledCentreline.GetPoint(line.GetId(ptId + 1)), \
                                       t, closestPoint)

            if tmpDistance <= distance:
                pos = ptId
                distance = tmpDistance
                lineId = tmpLineId

        line.Reset()
        tmpLineId = tmpLineId + 1

    # Change the sign of the segment number if it belongs to the parent branch.
    if lineId == 0:
        pos = pos - (lineLengths[lineId] - 1)

    return lineId, pos, t.get(), math.sqrt(distance)

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
    print centreline.GetNumberOfPoints(), 'points in the centreline...'

    # Resample the centreline to EC resolution.
    lineArray = centreline.GetLines()
    line = vtk.vtkIdList()
    points = vtk.vtkPoints()

    newCentrelinePoints = vtk.vtkPoints()
    newCentrelineLines = vtk.vtkCellArray()

    ecResolution = numAxialQuads * numECsPerCol

    # Iterate over all lines in the centreline.
    lineArray.InitTraversal()
    while lineArray.GetNextCell(line):

        # Get points just for this line.
        centreline.GetPoints().GetPoints(line, points)

        # Construct a parametric spline with these points.
        parametricSpline = vtk.vtkParametricSpline()
        parametricSpline.SetPoints(points)
        parametricSpline.SetLeftConstraint(2)
        parametricSpline.SetRightConstraint(2)

        # Construct a parametric line source based onthe spline.
        splinePointSource = vtk.vtkParametricFunctionSource()
        splinePointSource.SetParametricFunction(parametricSpline)
        splinePointSource.SetUResolution(ecResolution)
        splinePointSource.Update()

        tmpPoints = splinePointSource.GetOutput()
        newPolyLine = vtk.vtkPolyLine()

        # Insert the points from the resampled segment one-by-one into the new points
        # container, while updating the new line with the ids of each newly inserted point.
        for ptId in range(tmpPoints.GetNumberOfPoints()):
            newPolyLine.GetPointIds().InsertNextId(newCentrelinePoints.InsertNextPoint(tmpPoints.GetPoint(ptId)))

        newCentrelineLines.InsertNextCell(newPolyLine)

    # Create a new resampled centreline.
    resampledCentreline = vtk.vtkPolyData()
    resampledCentreline.SetPoints(newCentrelinePoints)
    resampledCentreline.SetLines(newCentrelineLines)

    # Transform to the nm scale space.
    transform = vtk.vtkTransform()
    transform.Scale((1000,1000,1000))

    transformFilter = vtk.vtkTransformFilter()
    transformFilter.SetTransform(transform)
    transformFilter.SetInput(resampledCentreline)
    transformFilter.Update()

    resampledCentreline = transformFilter.GetOutput()

    print resampledCentreline.GetNumberOfPoints(), 'points in the resampled centreline...'

    # Put the ecMesh through centroids filter.
    centroidFilter = vtk.vtkCellCenters()
    centroidFilter.SetInput(ecMesh)
    centroidFilter.Update()

    centroids = centroidFilter.GetOutput()

    # Iterate over each centroid and find the closest segment
    centroidPoints = centroids.GetPoints()

    atpArray = vtk.vtkDoubleArray()
    atpArray.SetName("initialATP")

    # Only for DEBUG output.
    tArray = vtk.vtkDoubleArray()
    tArray.SetName("T")

    # Only for DEBUG output.
    distArray = vtk.vtkDoubleArray()
    distArray.SetName("Dist")

    # Only for DEBUG output.
    branchArray = vtk.vtkDoubleArray()
    branchArray.SetName("Branch")

    # Only for DEBUG output.
    posArray = vtk.vtkDoubleArray()
    posArray.SetName("Pos")

    totalPoints = centroidPoints.GetNumberOfPoints()
    complete = 0
    for ptId in range(totalPoints):

        # Progress reporting.
        tmp = (ptId * 100) / float(totalPoints)
        if tmp - complete >= 1:
            complete = int(tmp)
            print complete, '% complete, processing point', ptId, 'out of', totalPoints, '...'
        branchId, pos, t, distance = getParametricDistance(resampledCentreline, centroidPoints.GetPoint(ptId))

        tArray.InsertNextValue(t)
        distArray.InsertNextValue(distance)
        branchArray.InsertNextValue(branchId)
        posArray.InsertNextValue(pos)

        atpArray.InsertNextValue(sigmoidATP(pos))

    print '100 % complete.'

    debugAtpDataset = ecMesh

    # Assert the number of cells is equal to the number of items in the cell arrays.
    assert atpArray.GetNumberOfTuples() == debugAtpDataset.GetNumberOfCells(), "Number of cells (%d) and cell data values (%d) mismatch." % (atpArray.GetNumberOfTuples(), debugAtpDataset.GetNumberOfCells())

    debugAtpDataset.GetCellData().AddArray(atpArray)
    debugAtpDataset.GetCellData().AddArray(tArray)
    debugAtpDataset.GetCellData().AddArray(distArray)
    debugAtpDataset.GetCellData().AddArray(branchArray)
    debugAtpDataset.GetCellData().AddArray(posArray)

    debugAtpMapWriter = vtk.vtkXMLPolyDataWriter()
    debugAtpMapWriter.SetFileName(debugAtpFile)
    debugAtpMapWriter.SetInput(debugAtpDataset)
    debugAtpMapWriter.Update()

    pointsToVerticesFilter = vtk.vtkVertexGlyphFilter()
    pointsToVerticesFilter.SetInput(centroids)
    pointsToVerticesFilter.Update()

    atpDataset = pointsToVerticesFilter.GetOutput()
    atpDataset.GetCellData().AddArray(atpArray)

    atpMapWriter = vtk.vtkXMLPolyDataWriter()
    atpMapWriter.SetFileName(atpFile)
    atpMapWriter.SetInput(atpDataset)
    atpMapWriter.Update()
    
    # Provide a quick visualisation of the ATP profile for validation.
    pointsX = range(-int(ecResolution - 1), int(ecResolution))
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