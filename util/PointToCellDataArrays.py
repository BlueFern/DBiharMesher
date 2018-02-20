'''
This script looks through a directory of .vtp files and converts the point
data to cell data for them all. The output is in a new directory.
'''

import glob
import os
import vtk

# This should really expect a directory as a parameter.

def pointToCellData():
    
    atpFiles = sorted(glob.glob('*.vtp'))
    
    if not atpFiles:
        exit("No atp files found")

    for inputFile in atpFiles: 

        print('Reading', inputFile)
        atpReader = vtk.vtkXMLPolyDataReader()
        atpReader.SetFileName(inputFile)
        atpReader.Update()

        atpDataset = atpReader.GetOutput()
        
        pointToCell = vtk.vtkPointDataToCellData()
        pointToCell.SetInputData(atpDataset)
        pointToCell.PassPointDataOn()
        pointToCell.Update()
        
        convertedData = pointToCell.GetOutput()
        
        # Remove the point data arrays, they exist as cell data.
        convertedData.GetPointData().RemoveArray('ATP')
        convertedData.GetPointData().RemoveArray('tau')
        convertedData.GetCellData().RemoveArray('p')
    
        atpWriter = vtk.vtkXMLPolyDataWriter()
        atpWriter.SetInputData(convertedData)
        atpWriter.SetFileName("_" + inputFile)
        atpWriter.Update()


def usage():
    print('This script is to be run with global parameters (angle) set in the calling script.')

if __name__ == '__main__':
    print("Starting", os.path.basename(__file__))
    usage()
    print("Exiting", os.path.basename(__file__))
