# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 14:38:16 2015

Create an initial JPLC Profile for an ECs Mesh.
"""

import os
import sys
import vtk

numBranches = 3
numQuads = 72
numECsPerCol = 4

meshFile = "quadMeshFullECc216.vtp"

# TODO: To make this script to be able to process meshes of arbitrary tree
# complexity, the script needs the following:
# 1. Read command line parameter for the mesh.
# 2. Mesh needs to have its branches labelled with unique itentifiers.

def main():
    # Report our CWD just for testing purposes.
    print "Current working directory:", os.getcwd()    

    # Read in the mesh.
    meshReader = vtk.vtkXMLPolyDataReader()
    meshReader.SetFileName(meshFile)
    meshReader.Update()

    mesh = meshReader.GetOutput()
    
    gridCoords = mesh.GetCellData().GetArray("gridCoords")
    
    # This is the part that will have to change, when the branches are labelled.
    
    
    # Iterate over all quad cells in the mesh and pull out the first component of the
    # "gridCoords" cell data array.
    # Here we assuming the mesh contains only quad cells.
    for cellId in range(0, mesh.GetNumberOfCells()):
        coords = gridCoords.GetTuple(cellId)
        print coords[0]
    

    # Put it through centroids filter.
    # Use VTK centroid filter to get the centroids in the right order
    # from the reorderedECMeshBranch.
    centroidFilter = vtk.vtkCellCenters()
    centroidFilter.SetInput(mesh)
    # centroidFilter.Update()

    # Create a vertex cell for each point.
    pointsToVerticesFilter = vtk.vtkVertexGlyphFilter()
    pointsToVerticesFilter.SetInputConnection(centroidFilter.GetOutputPort())
    #pointsToVerticesFilter.SetInput(centroidFilter.GetOutput())
    pointsToVerticesFilter.Update()
    
    centroidsAndVertices = pointsToVerticesFilter.GetOutput()
    

    
    
    
    verticesArray = centroidsAndVertices.GetVerts()
    verticesArray.InitTraversal()
    

# Create a new cell data array "parametricDistance" in the centroids dataset
# with the values from the previous step. This array is really not needed for
# any computations is good to have for verification. Also, remember the max
# value for each branch.

# Iterate over the "parametricDistance" cell array and for each value call the
# agonistValue function and store the computed values in the new "InitialJPLC"
# cell data array.    
    


if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)