#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# This scipt grabs a cell or point array from the first input dataset, which in our case would
# be an ATP map where the data is associated with cell centre vertices (cell data), or the points
# at those vertices (point data), and converts them to cell data associated with cells, which
# in our case represend EC cells.

import os
import sys
import vtk
import argparse

def CopyData():
    print args
    input_1 = vtk.vtkPolyData()

    # Read input_1.
    if os.path.splitext(args.input_1)[1] == ".vtp":
        reader1 = vtk.vtkXMLPolyDataReader()
        reader1.SetFileName(args.input_1)
        reader1.Update()
        input_1 = reader1.GetOutput()
    else:
        raise Exception("Expected vtp as input_1.")

    # print input_1

    # Read input_2
    if os.path.splitext(args.input_2)[1] == ".vtp":
        reader2 = vtk.vtkXMLPolyDataReader()
        reader2.SetFileName(args.input_2)
        reader2.Update()
        input_2 = reader2.GetOutput()
    else:
        raise Exception("Expected vtp as input_2.")

    # Names for arrays to copy.
    point_arrays = []
    cell_arrays = []

    for array_name in args.array_names:
        if input_1.GetPointData().HasArray(array_name):
            point_arrays.append(array_name)
        else:
            print array_name, "not found in point data arrays"

        if input_1.GetCellData().HasArray(array_name):
            cell_arrays.append(array_name)
        else:
            print array_name, "not found in cell data arrays"

    # These are the names of point arrays to copy.
    print "Point data arrys to copy", point_arrays
    print "Cell data arrays to copy", cell_arrays

    # Copy point arrays and save the output.
    output = vtk.vtkPolyData()
    output.DeepCopy(input_2)

    # Copy point data arrays.
    for array_name in point_arrays:
        array = input_1.GetPointData().GetArray(array_name)

        # Make sure the number of cells in the output matches the number of tuples in the source array.
        assert (array.GetNumberOfTuples() == output.GetNumberOfCells()), "Number of tuples in array from input_1 does not match the number of cells in input_2."
        output.GetCellData().AddArray(array)

    # Copy cell data arrays.
    for array_name in cell_arrays:
        array = input_1.GetCellData().GetArray(array_name)

        # Make sure the number of cells in the output matches the number of tuples in the source array.
        assert (array.GetNumberOfTuples() == output.GetNumberOfCells()), "Number of tuples in array from intput_1 does not match the number of cells in input_2."
        output.GetCellData().AddArray(array)

    # Write output.
    if os.path.splitext(args.output)[1] == ".vtp":
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(args.output)
        writer.SetInput(output)
        writer.Update()
    else:
        raise Exception("Expected vtp as output.")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Grab point or cell data array from first input and copy them into second input.")
    parser.add_argument("input_1", help="First input filename")
    parser.add_argument("input_2", help="Second input filename")
    parser.add_argument("output", help="Output filename")
    parser.add_argument("array_names", nargs="+", help="Names of arrays to copy")
    args = parser.parse_args()


    CopyData()




