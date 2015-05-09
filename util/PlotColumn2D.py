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

suffix = ''
yLine = None

def plotColumn2D(fileList, splitPoint = 0, startPoint = 0, endPoint = -1):
    # Report our CWD just for testing purposes.
    print "CWD:", os.getcwd()    

    sortNicely(fileList)

    allRows = []

    for file in fileList:
        print 'Reading', file

        figName = os.path.split(os.getcwd())[1]
        if suffix != '':
            figName = figName + '.' + suffix

        file_name = os.path.abspath(file)
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(file_name)
        reader.Update()
        data = reader.GetOutput().GetCellData().GetArray(0)

        row = []
        for i in range(data.GetNumberOfTuples()):
            row.append(data.GetValue(i))

        # Fix the croocked backwards ordering in the output data. Sigh...
        # Also, splice the branches.
        if isinstance(splitPoint, str):
            if splitPoint == 'mid':
                splitPoint = len(row) / 2
            else:
                print 'Can\'t split row with splitPoint =', splitPoint
                return

        r1 = row[:splitPoint]
        r2 = row[len(r1):]
        row = r2 + r1

        # Trim.
        row = row[startPoint:endPoint]
        row.reverse()

        allRows.append(row)

    array2D = numpy.array(allRows)
    array2D = numpy.transpose(array2D)

    plt.pcolormesh(array2D)

    if yLine != None:
        plt.axhline(yLine, color='r')
    
    plt.title(figName)
    plt.xlabel('Time (sec.)')
    plt.ylabel('Cell (ord.)')
    plt.colorbar()
    plt.axis([0, array2D.shape[1], 0, array2D.shape[0]])
    plt.tight_layout()

    plt.savefig(figName + '.png', bbox_inches='tight', dpi=400)
    plt.show()

def usage():
    print "This script is to be run with global parameters (input files, splitPoint, startPoint, endPoint) set in the calling script."

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    usage()
    print "Exiting", os.path.basename(__file__)
