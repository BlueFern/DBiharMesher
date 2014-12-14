#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# SynthCentrelineGenerator.py by Constantine Zakkaroff
#
# This scipt builds a centreline tree from a nested list representing the tree topology.
# The fist item in each list is the domain, followed by the left
# and right branches if any.
#
# Output centrelines are saved in VTK legacy format in the current working directory
# in a file "centreline.vtk".
#
# Parameters to control the generated tree points are to be set in this file.
# VTK visualisation of the tree is shown at the end.
# 
# Tree grows from (0,0,0) towards positive x-axis. Bifurcations occur at 45 degrees
# by defualt, or at a given angle if specfied.
#
# Domains are either given as a number, which represents the domain length, or a 
# tuple that represents (domain length, angle).
#
# All angles are measured clockwise from the positive y-axis, and are specified
# in degrees.
# Modified by: Stewart Dowding 15/12/14

# TODO:
# Name and location of the outputfile should be specified from keyboard input at runtime.

import sys
import vtk
import math

# Lists to be used as input. A specific list is passed as the argument to
# buildCentreline in main.

branchAngle = math.pi / 4.0

segmentList0 = [20,[(20,0),None,None],[(20,45),None,None]]
segmentList1 = [16,[8,[4,[2,None,None],[2,None,None]],None],[8,[4,None,None],None]]
segmentList2 = [20,[(20,60),[(20,80),None,None],[20,None,None]],[(20,135),None,None]]

# Paremeters specifying the properties of the generated centrelines.

scaling = 1.0
step = 0.1
radiusBase = 2.0
radiusDelta = 0.0025

# If sphereRadius is set to None, the centreline is generated in the XY plane.
# Otherwise the centreline is wrapped on a sphere of the specified radius.
sphereRadius = None

points = vtk.vtkPoints()
lines = vtk.vtkCellArray()
radii = vtk.vtkDoubleArray()

def buildCentreline(segmentList, firstId = 0, firstPt = (0.0,0.0,0.0), direction = 0.0, radius = radiusBase):
    print 'Processing centreline:', segmentList
    
    domain = None

    try:
        domain = segmentList[0]
    except IndexError:
        sys.exit("Missing domain length in list " + str(list) + ".")
        
    if isinstance(domain, tuple):
        domainLength = domain[0]
        angle = math.radians(domain[1])
    else:    
        if not isinstance(domain, (int, float, long)):
            sys.exit("Domain length should be a number, not a(an) " + str(type(domainLength)) + ".")
        angle = branchAngle
        domainLength = domain
    
    leftBranch = None
    try:
        leftBranch = segmentList[1]
    except IndexError:
        pass
    
    rightBranch = None
    try:
        rightBranch = segmentList[2]
    except IndexError:
        pass
    
    line = vtk.vtkPolyLine()
    
    if  direction < 0:
        angle = math.pi - angle
    
    for pId in range (0, int(domainLength / step) + 1):
            nextPt = ((firstPt[0] + (1.0 if direction == 0.0 else math.sin(angle)) * (step * pId)) * scaling, \
                        (firstPt[1] + (1.0 if direction == 0.0 else math.cos(angle)) * (step * pId) * direction) * scaling , \
                        firstPt[2] * scaling)
                        
            if sphereRadius != None:
                longitude = nextPt[0] / sphereRadius
                latitude = 2 * math.atan(math.exp(nextPt[1] / sphereRadius)) - math.pi / 2.0
                nextPtS = (sphereRadius * math.cos(latitude) * math.cos(longitude),\
                            sphereRadius * math.cos(latitude) * math.sin(longitude),\
                            sphereRadius * math.sin(latitude))
            
                if firstId == 0 or pId > 0:
                    nextId = points.InsertNextPoint(nextPtS)
                else:
                    nextId = firstId
            else:
                if firstId == 0 or pId > 0:
                    nextId = points.InsertNextPoint(nextPt)
                else:
                    nextId = firstId
            
            line.GetPointIds().InsertNextId(nextId)
            
            if firstId == 0 or pId > 0:
                radiusVal = radius - (radiusDelta * pId)
                radii.InsertNextTuple((radiusVal,))
    
    lines.InsertNextCell(line)
    
    if isinstance(leftBranch, list):
        buildCentreline(leftBranch, nextId, nextPt, 1.0, radiusVal)

    if isinstance(rightBranch, list):
        buildCentreline(rightBranch, nextId, nextPt, -1.0, radiusVal)

def main():
    buildCentreline(segmentList2)
    
    print "Number of points in the centreline:", points.GetNumberOfPoints()
    
    centreline = vtk.vtkPolyData()
    centreline.SetPoints(points)
    centreline.SetLines(lines)
    centreline.GetPointData().SetScalars(radii)


    if sphereRadius != None:
        origin = centreline.GetPoint(0)
        transform = vtk.vtkTransform()
        transform.Translate(-origin[0],-origin[1],-origin[2])
        transformFilter = vtk.vtkTransformPolyDataFilter()
        transformFilter.SetInput(centreline)
        transformFilter.SetTransform(transform)
        
        transformFilter.Update()
        centreline = transformFilter.GetOutput()
    
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
    style = vtk.vtkInteractorStyleTrackballCamera()
    interactor.SetInteractorStyle(style)
    interactor.SetRenderWindow(rendererWindow)
        
    interactor.Initialize()
    interactor.Start()
    
    rendererWindow.Finalize()
    interactor.TerminateApp()

if __name__ == '__main__':
    main()