# -*- coding: utf-8 -*-
"""
Use a 'reference' ATP/Jplc map for a symmetric c4080 bifurcation to be applied to a given EC mesh to generate
an ATP map specific to some geometry for an asymmetric c4080 bifurcation.
"""

import os
import vtk

# Run in the directory where this script is.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def copyATPFromRef(angle):
    print "Processing angle", angle, "..."
    # Read reference file.
    # Get Jplc double array from point data.
    referenceReader = vtk.vtkXMLPolyDataReader()
    referenceReader.SetFileName("quadMeshFullATPc4080_reference.vtp")
    referenceReader.Update()
    
    Jplc = referenceReader.GetOutput().GetCellData().GetArray("initialATP")
    
    # Read EC mesh.
    # Get centroids.
    ecMeshReader = vtk.vtkXMLPolyDataReader()
    ecMeshReader.SetFileName(os.path.join(str(angle), "quadMeshFullECc4080.vtp"))
    ecMeshReader.Update()
    
    ecSurface = ecMeshReader.GetOutput()
    
    ecCentroidFilter = vtk.vtkCellCenters()
    ecCentroidFilter.VertexCellsOn()
    ecCentroidFilter.SetInput(ecSurface)
    ecCentroidFilter.Update()
    
    ecCentroids = ecCentroidFilter.GetOutput()

    # Write EC mesh with Jplc as cell data.
    ecSurface.GetCellData().AddArray(Jplc)
    atpSurfaceWriter = vtk.vtkXMLPolyDataWriter()
    atpSurfaceWriter.SetInput(ecSurface)
    atpSurfaceWriter.SetFileName(os.path.join(str(angle), "quadMeshFullATPc4080_.vtp"))
    atpSurfaceWriter.Update()

    # Write EC centroids with Jpls as cell data.
    ecCentroids.GetCellData().AddArray(Jplc)
    atpPointsWriter = vtk.vtkXMLPolyDataWriter()
    atpPointsWriter.SetInput(ecCentroids)
    atpPointsWriter.SetFileName(os.path.join(str(angle), "quadMeshFullATPc4080.vtp"))
    print atpPointsWriter.GetFileName()
    atpPointsWriter.Update()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    # Set angles and directories and call copyATPFromRef.
    for angle in range(60, 120, 10):
        copyATPFromRef(angle)
    print "Exiting", os.path.basename(__file__)