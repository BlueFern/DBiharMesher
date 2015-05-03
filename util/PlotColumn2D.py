# -*- coding: utf-8 -*-
"""
Read SMC Ca2+ values from a line of cells in temporal series and produce a 2D plot
"""
import re
import os
import vtk
import numpy
import matplotlib.pyplot as plt

def tryInt(s):
    try:
        return int(s)
    except:
        return s

def alphaNumKey(s):
    """ Turn a string into a list of string and number chunks.
        "z23a" -> ["z", 23, "a"]
    """
    return [tryInt(c) for c in re.split('([0-9]+)', s)]

def sortNicely(l):
    """ Sort the given list in the way that humans expect.
    """
    l.sort(key=alphaNumKey)

# This is to be done in the calling script.
# import glob
# file_list = glob.glob("*.vtu")

def plotColumn2D(fileList, splitPoint, startPoint, endPoint):
    sortNicely(fileList)
  
    allRows = []
    
    for file in fileList:
        file_name = os.path.abspath(file)
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(file_name)
        reader.Update()
        data = reader.GetOutput().GetCellData().GetArray(0)
    
        row = []
        for i in range(data.GetNumberOfTuples()):
            row.append(data.GetValue(i))
        
        # Fixe croocked backwards ordering in the output data. Sigh...
        r1 = row[:splitPoint]
        r2 = row[splitPoint + 1:]
        row = r2 + r1
        
        # Trim.
        row = row[startPoint:endPoint]
        row.reverse()
    
        allRows.append(row)
    
    array2D = numpy.array(allRows)
    array2D = numpy.transpose(array2D)
    
    plt.pcolor(array2D)
    plt.show()

def main():
    print "This script is to be run with global parameters (input files, splitPoint, startPoint, endPoint) set in the calling script."

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)