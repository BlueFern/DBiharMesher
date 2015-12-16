import vtk
import math
import numpy as np

atp_min = 0.2
atp_max = 0.7

if __name__ == "__main__":

    # Read ECs.
    ecReader = vtk.vtkXMLPolyDataReader()
    ecReader.SetFileName('quadMeshFullECc256.vtp')
    ecReader.Update()

    # Rememeber the cells.
    polydata = ecReader.GetOutput()

    # Create new ATP map array.
    newMapArray = vtk.vtkDoubleArray()
    newMapArray.SetName('initialATP')
    newMapArray.SetNumberOfComponents(1)
    newMapArray.SetNumberOfTuples(polydata.GetNumberOfCells())

    # Get cell centres.
    centresFilter = vtk.vtkCellCenters()
    centresFilter.SetInput(polydata)
    centresFilter.Update()

    centres = centresFilter.GetOutput()

    _, _, y_min, y_max, _, _ = centres.GetBounds()
    y_range = y_max - y_min
    atp_range = atp_max - atp_min

    # Set background value in the ATP array to vary linearly from atp_min to atp_max.
    # Create rectangular area of high ATP in the middle.
    for cid in range (0, centres.GetNumberOfPoints()):
        x, y, _ = centres.GetPoint(cid)

        # Rescale to range.
        atp_val = (y - y_min) * (atp_range) / (y_range) + atp_min

        newMapArray.SetValue(cid, atp_val)

    # Set new map.
    polydata.GetCellData().AddArray(newMapArray)
    polydata.GetCellData().SetActiveScalars('initialATP')

    # Write the map.
    mapWriter = vtk.vtkXMLPolyDataWriter()
    mapWriter.SetFileName('quadMeshFullATPc256.vtp')
    mapWriter.SetInput(polydata)
    mapWriter.Update()

    # Setup actor and mapper
    mapper = vtk.vtkPolyDataMapper()
    if vtk.VTK_MAJOR_VERSION <= 5:
        mapper.SetInput(polydata)
    else:
        mapper.SetInputData(polydata)
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    # Setup render window, renderer, and interactor
    renderer = vtk.vtkRenderer()
    renderWindow = vtk.vtkRenderWindow()
    renderWindow.AddRenderer(renderer)
    renderWindowInteractor = vtk.vtkRenderWindowInteractor()
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderer.AddActor(actor)
    renderWindow.Render()
    renderWindowInteractor.Start()


