'''
This script looks through a directory of .vtp files and converts the point
data to cell data for them all. The output is in a new directory.
'''

import glob
import os
import vtk

angle = 0

def pointToCellData():
    
    print "Using angle: " + str(angle)
    atpFiles = glob.glob(str(angle) + '/*.vtp')
    if not atpFiles:
        exit("No atp files found")
        
    if not os.path.exists(str(angle) + '_new'):
        os.makedirs(str(angle) + '_new')
    
    for inputFile in atpFiles: 

        print 'Reading', inputFile
        atpReader = vtk.vtkXMLPolyDataReader()
        atpReader.SetFileName(inputFile)
        atpReader.Update()
    
        atpDataset = atpReader.GetOutput()
        
        pointToCell = vtk.vtkPointDataToCellData()
        pointToCell.SetInputData(atpDataset)
        pointToCell.PassPointDataOn()
        pointToCell.Update()
        
        newArray = pointToCell.GetOutput()
        
        # Remove the point data arrays, they exist as cell data
        newArray.GetPointData().RemoveArray('ATP')
        newArray.GetPointData().RemoveArray('tau')
        newArray.GetCellData().RemoveArray('p')
    
        atpWriter = vtk.vtkXMLPolyDataWriter()
        atpWriter.SetInputData(newArray)
        atpWriter.SetFileName(str(angle) + '_new/' + inputFile[3:])
        atpWriter.Update()


def usage():
    print 'This script is to be run with global parameters (input ATP, angle) set in the calling script.'

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    usage()
    print "Exiting", os.path.basename(__file__)