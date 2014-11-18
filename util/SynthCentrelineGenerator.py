#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# by Constantine Zakkaroff
#
# This scipt builds a tree from a nested list representing the tree topology.
# The fist item in each list is the length of the branch, followed by the left
# and right branches if any.
#
# Parameters to control the generated tree points are to be set in this file.
# VTK visualisation of the tree is shown at the end.

import sys
import vtk
import math

treeList0 = [20,[20,None,None],[20,None,None]]
treeList1 = [16,[8,[4,[2,None,None],[2,None,None]],None],[8,[4,None,None],None]]

branchAngle = math.pi / 4.0
scaling = 1.0
step = 0.1
radiusBase = 2.0

points = vtk.vtkPoints()
lines = vtk.vtkCellArray()
radii = vtk.vtkDoubleArray()

def buildTree(treeList, startingPoint = (0.0,0.0,0.0), direction = 0.0, radius = radiusBase):
    print treeList
    
    domainLength = None
    try:
        domainLength = treeList[0]
    except IndexError:
        sys.exit("Missing domain length in list " + str(list) + ".")
    if not isinstance(domainLength, (int, float, long)):
        sys.exit("Domain length should be a number, not a(an) " + str(type(domainLength)) + ".")
 
    leftBranch = None
    try:
        leftBranch = treeList[1]
    except IndexError:
        pass

    rightBranch = None
    try:
        rightBranch = treeList[2]
    except IndexError:
        pass

    line = vtk.vtkPolyLine()

    firstPoint = None
    firstPointId = None
  
    if isinstance(startingPoint, tuple):
        firstPoint = startingPoint
        firstPointId = points.InsertNextPoint(firstPoint)
        radii.InsertNextTuple((radius,))
    elif isinstance(startingPoint, (int, long)):
        firstPointId = startingPoint
        firstPoint = points.GetPoint(firstPointId)
    else:
        sys.exit("Starting point can be a tuple or an integer only, got " + str(type(startingPoint)) + ".")

    line.GetPointIds().InsertNextId(firstPointId)

    nextPointId = None
    for p in range (1, int(domainLength / step) + 1):
            nextPoint = ((firstPoint[0] + (1.0 if direction == 0.0 else math.sin(branchAngle)) * (step * p)) * scaling, \
                        (firstPoint[1] + (1.0 if direction == 0.0 else math.cos(branchAngle)) * (step * p) * direction) * scaling, \
                        firstPoint[2] * scaling)
            nextPointId = points.InsertNextPoint(nextPoint)
            line.GetPointIds().InsertNextId(nextPointId)
            print 'p: (%.3f, %.3f, %.3f)' % nextPoint,
            radius = radius - (step / 40.0)
            print 'r: %.3f' % radius
            radii.InsertNextTuple((radius,))
    
    lines.InsertNextCell(line)

    if leftBranch != None:
        buildTree(leftBranch, nextPointId, -1.0, radius)
            
    if rightBranch != None:
        buildTree(rightBranch, nextPointId, 1.0, radius)

def main():
    print "Starting..."

    buildTree(treeList0)
    
    print "Number of points in the tree:", points.GetNumberOfPoints()
    
    centreline = vtk.vtkPolyData()
    centreline.SetPoints(points)
    centreline.SetLines(lines)
    centreline.GetPointData().SetScalars(radii)
    
    writer = vtk.vtkPolyDataWriter()
    writer.SetInput(centreline)
    writer.SetFileName("centreline.vtk")
    writer.SetFileTypeToASCII()
    writer.Write()
    
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInput(centreline)
    
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    
    renderer = vtk.vtkRenderer()
    renderer.AddActor(actor)
    
    rendererWindow = vtk.vtkRenderWindow()
    rendererWindow.AddRenderer(renderer)
    
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(rendererWindow)
    interactor.Initialize()
    interactor.Start()
    
    rendererWindow.Finalize()
    interactor.TerminateApp()
    
    print "Done."

if __name__ == '__main__':
    main()