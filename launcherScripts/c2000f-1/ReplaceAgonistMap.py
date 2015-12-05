import vtk

bgATp = 0.35

ecReader = vtk.vtkXMLPolyDataReader()
ecReader.SetFileName('quadMeshFullECc2000.vtp')
ecReader.Update()

polydata = ecReader.GetOutput()

# Create new ATP map array.
newMapArray = vtk.vtkDoubleArray()
newMapArray.SetName('initialATP')
newMapArray.SetNumberOfComponents(1)
newMapArray.SetNumberOfTuples(polydata.GetNumberOfCells())

# Set background value in the ATP array.
for cid in range(0, polydata.GetNumberOfCells()):
    newMapArray.SetValue(cid, bgATp)

# Get cell centres.
centresFilter = vtk.vtkCellCenters()
centresFilter.SetInput(polydata)
centresFilter.Update()

centres = centresFilter.GetOutput()

rH = 1000
pointH = (5000, 7507.5, 0)

# Do HI ATP by running over all cell centres in polydata and processing only the ones that fall in the specified radius.
for cid in range (0, centres.GetNumberOfPoints()):
    if vtk.vtkMath.Distance2BetweenPoints(pointH, centres.GetPoint(cid)) < rH**2:
        newMapArray.SetValue(cid, 0.68)

rL1 = 1000
pointL1 = (3000, 3022.5, 0)

# Do LO ATP 1 by running over all cell centress in polydata and processing only the ones that fall in the specified radius.
for cid in range (0, centres.GetNumberOfPoints()):
    if vtk.vtkMath.Distance2BetweenPoints(pointL1, centres.GetPoint(cid)) < rL1**2:
        newMapArray.SetValue(cid, 0.25)

rL2 = 1000
pointL2 = (7000, 3022.5, 0)

# Do LO ATP 2 by running over all cell centres in polydata and processing only the ones that fall in the specified radius.
for cid in range (0, centres.GetNumberOfPoints()):
    if vtk.vtkMath.Distance2BetweenPoints(pointL2, centres.GetPoint(cid)) < rL2**2:
        newMapArray.SetValue(cid, 0.3)

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

