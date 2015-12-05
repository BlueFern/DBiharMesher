# -*- coding: utf-8 -*-
"""
Created ATP profile input data based on rescaling CDF output.

"""

import os
import vtk

inputFile = ''
inputSurfaceFile = ''
outputFile = ''
outputSurfaceFile = ''
outMin = 0.0
outMax = 1.0

def rescale(val, inMin, inMax):
    return (val - inMin) * (outMax - outMin) / (inMax - inMin) + outMin

def rescaleATPToRange():
    # This is where the data is for testing purposes.
    print "Current working directory:", os.getcwd()

    print 'Reading', inputFile
    atpReader = vtk.vtkXMLPolyDataReader()
    atpReader.SetFileName(inputFile)
    atpReader.Update()

    atpDataset = atpReader.GetOutput()
    atp = atpDataset.GetPointData().GetArray('ATP')

    inMin, inMax = atp.GetRange()

    atpR = vtk.vtkFloatArray()
    atpR.SetName('initialATP')    
    
    for i in range(atp.GetNumberOfTuples()):
        val = atp.GetTuple(i)[0]
        atpR.InsertNextValue(rescale(val, inMin, inMax))

    atpDataset.GetCellData().AddArray(atpR)
    
    atpDataset.GetPointData().RemoveArray('ATP')
    atpDataset.GetPointData().RemoveArray('tau_w')
    atpDataset.GetCellData().RemoveArray('initialJPLC')

    print 'Writing', outputFile
    atpWriter = vtk.vtkXMLPolyDataWriter()
    atpWriter.SetInput(atpDataset)
    atpWriter.SetFileName(outputFile)
    atpWriter.Update()

    print 'Reading', inputSurfaceFile
    surfaceReader = vtk.vtkXMLPolyDataReader()
    surfaceReader.SetFileName(inputSurfaceFile)
    surfaceReader.Update()

    surface = surfaceReader.GetOutput()
    surface.GetCellData().AddArray(atpR)

    print 'Writing', outputSurfaceFile
    surfaceWriter = vtk.vtkXMLPolyDataWriter()
    surfaceWriter.SetInput(surface)
    surfaceWriter.SetFileName(outputSurfaceFile)
    surfaceWriter.Update()

def usage():
    print 'This script is to be run with global parameters (input ATP, desired value range, etc.) set in the calling script.'

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    usage()
    print "Exiting", os.path.basename(__file__)