# -*- coding: utf-8 -*-
"""
Created on Fri May  1 12:09:31 2015

@author: constantine
"""

import os
import vtk

inputFile = 'quadMeshFullJPLCc4080.vtp'
outputFile = 'quadMeshFullATPc4080.vtp'

outMin = 0.1
outMax = 1.1

def rescale(val, inMin, inMax):
    return (val - inMin) * (outMax - outMin) / (inMax - inMin) + outMin

def main():
    # This is where the data is for testing purposes.
    print "Current working directory:", os.getcwd()

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
        
    atpWriter = vtk.vtkXMLPolyDataWriter()
    atpWriter.SetInput(atpDataset)
    atpWriter.SetFileName(outputFile)
    atpWriter.Update()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)