# -*- coding: utf-8 -*-
"""
Generate centreline and write it out as .vtk legacy format.
"""

import os
import vtk

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

def main():
    # Read reference file.
    # Get Jplc double array from point data.
    referenceReader = vtk.vtkXMLPolyDataReader()
    referenceReader.SetFileName("quadMeshFullATPc4080_reference.vtp")
    referenceReader.Update()
    
    Jplc = referenceReader.GetOutput().GetCellData().GetArray("initialATP")
    
    # Read EC mesh.
    # Get centroids.
    ecMeshReader = vtk.vtkXMLPolyDataReader()
    ecMeshReader.SetFileName("quadMeshFullECc4080.vtp")
    ecMeshReader.Update()
    
    ecSurface = ecMeshReader.GetOutput()
    
    ecCentroidFilter = vtk.vtkCellCenters()
    ecCentroidFilter.VertexCellsOn()
    ecCentroidFilter.SetInput(ecSurface)
    ecCentroidFilter.Update()
    
    ecCentroids = ecCentroidFilter.GetOutput()

    # Write EC mesh with Jplc as cell data.
    # Write EC centroids with Jpls as cell data.
    ecSurface.GetCellData().AddArray(Jplc)
    atpSurfaceWriter = vtk.vtkXMLPolyDataWriter()
    atpSurfaceWriter.SetInput(ecSurface)
    atpSurfaceWriter.SetFileName("quadMeshFullATPc4080_.vtp")
    atpSurfaceWriter.Update()
    
    ecCentroids.GetCellData().AddArray(Jplc)
    atpPointsWriter = vtk.vtkXMLPolyDataWriter()
    atpPointsWriter.SetInput(ecCentroids)
    atpPointsWriter.SetFileName("quadMeshFullATPc4080.vtp")
    atpPointsWriter.Update()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)