import vtk
import math
import numpy as np

ecReader = vtk.vtkXMLPolyDataReader()
ecReader.SetFileName('quadMeshFullECc2000.vtp')
ecReader.Update()

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

atp_min = 0.2
atp_max = 0.7

# High ATP rectange coordinates.
r_x1 = 3000
r_x2 = 7000
r_y1 = 6000
r_y2 = 7000

_, _, y_min, y_max, _, _ = centres.GetBounds()
y_range = y_max - y_min

# Set background value in the ATP array to vary linearly from atp_min to atp_max.
# Create rectangular area of high ATP in the middle.
for cid in range (0, centres.GetNumberOfPoints()):
    x, y, _ = centres.GetPoint(cid)

    if (r_x1 < x and x < r_x2) and (r_y1 < y and y < r_y2):

        # ATP max value.
        atp_val = atp_max

    else:

        # Rescale to range.
        atp_val = (y - y_min) * (atp_max - atp_min) / (y_max - y_min) + atp_min

    newMapArray.SetValue(cid, atp_val)

# Set new map.
polydata.GetCellData().AddArray(newMapArray)
polydata.GetCellData().SetActiveScalars('initialATP')

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

mapWriter = vtk.vtkXMLPolyDataWriter()
mapWriter.SetFileName('quadMeshFullATPc2000.vtp')
mapWriter.SetInput(polydata)
mapWriter.Update()

