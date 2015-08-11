# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 12:39:48 2015

A horrible quick and dirty script to build a mesh in our desired format.

@author: sed59
"""
import vtk

xQuads = 4
yQuads = 2
quadLength = 200
quadHeight = 300

def buildATPMesh(polydata, filename):
    
    centroidFilter = vtk.vtkCellCenters()
    centroidFilter.VertexCellsOn()
    centroidFilter.SetInput(polydata)
    
    newPolydata = vtk.vtkPolyData()
    newPolydata = centroidFilter.GetOutput()
    centroidFilter.Update()
    
    ATPValues = vtk.vtkDoubleArray()
    ATPValues.SetName("initialATP")
    
    _, _, yMin, yMax, _, _ = polydata.GetBounds()
    yRange = yMax - yMin
    
    for pointId in range(0, newPolydata.GetNumberOfPoints()):
        _, y, _ = newPolydata.GetPoint(pointId)
        ATPValue = y / (yRange * 1.0)
        ATPValues.InsertNextValue(ATPValue)
        
    newPolydata.GetCellData().SetScalars(ATPValues)
    
    polyDataWriter = vtk.vtkXMLPolyDataWriter()
    polyDataWriter.SetFileName(filename)
    polyDataWriter.SetInput(newPolydata)    
    polyDataWriter.Write()


def buildMesh(xNumCells, yNumCells, filename):
    
    polydata = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    quads = vtk.vtkCellArray()

    counter = 0 
    
    xBase = 0
    yBase = 0
    
    xStep = quadLength/(xNumCells * 1.0) # Focing float division (python 2).
    yStep = quadHeight/(yNumCells * 1.0)
    
    branchId = vtk.vtkDoubleArray()
    branchId.SetName("branchId")

    
    # Three rectangles in space.
    for k in range(0, 3):
        if k == 1:
            xBase = - quadLength * (xQuads / 2)
            yBase = quadHeight * yQuads
            
        elif k == 2:
            xBase = + quadLength * (xQuads / 2) 
            yBase = quadHeight * yQuads     
                       
        # The framework for each quad. The total number is set by globals.
        for i in range(0, yQuads):
            for j in range(0, xQuads):
                
                    # For how ever many cells are within each main quad.
                    # The total number is decided by the function parameters.
                    for x in range(0, yNumCells):
                        for y in range(0, xNumCells):
                            
                            quad = vtk.vtkQuad()
                                                        
                            p0 = [y * xStep + xBase + (quadLength * j), 
                                  x * yStep + yBase + (quadHeight * i), 0]
                                  
                            p1 = [(y + 1) * xStep + xBase + (quadLength * j),
                                  x * yStep + yBase + (quadHeight * i), 0]
                                  
                            p2 = [(y + 1) * xStep + xBase + (quadLength * j), 
                                  (x + 1) * yStep + yBase + (quadHeight * i), 0]
                                  
                            p3 = [y * xStep + xBase + (quadLength * j),
                                  (x + 1) * yStep + yBase + (quadHeight * i), 0]
                                  
                            
                            points.InsertNextPoint(p0)
                            points.InsertNextPoint(p1)
                            points.InsertNextPoint(p2)
                            points.InsertNextPoint(p3)
                            
                            quad.GetPointIds().InsertId(0, counter)
                            quad.GetPointIds().InsertId(1, counter + 1)
                            quad.GetPointIds().InsertId(2, counter + 2)
                            quad.GetPointIds().InsertId(3, counter + 3)
                            counter += 4
                            
                            quads.InsertNextCell(quad)
                            branchId.InsertNextValue(k)

    polydata.SetPoints(points)
    polydata.SetPolys(quads)
    polydata.GetCellData().SetScalars(branchId)
    
    polyDataWriter = vtk.vtkXMLPolyDataWriter()
    polyDataWriter.SetFileName(filename)
    polyDataWriter.SetInput(polydata)    
    polyDataWriter.Write()

    return polydata    
    