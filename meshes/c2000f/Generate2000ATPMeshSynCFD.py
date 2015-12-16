import vtk
import math
import numpy as np

sigmoid_grad = -3.0

def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-sigmoid_grad * x))

inner_par_rad = 1/2.0
atp_base = 0.35

sigmoind_domain_min = -2.0
sigmoind_domain_max = 2.0
sigmoid_domain = sigmoind_domain_max - sigmoind_domain_min
scaling_to_sigmoid_domain = sigmoid_domain / (1.0 - inner_par_rad)

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
    newMapArray.SetValue(cid, atp_base)

# Get cell centres.
centresFilter = vtk.vtkCellCenters()
centresFilter.SetInput(polydata)
centresFilter.Update()

centres = centresFilter.GetOutput()

sphere_rad = 1500
centre = (5000, 6500, 0)
max_atp = 0.68

# Do HI ATP by running over all cell centres in polydata and processing only the ones that fall in the specified radius.
for cid in range (0, centres.GetNumberOfPoints()):

    # Distance from the point to sphere centre.
    dist2 = vtk.vtkMath.Distance2BetweenPoints(centres.GetPoint(cid), centre)

    # If the point falls into the sphere radius.
    if dist2 < sphere_rad**2:

        # Parametric distance from point to sphere centre.
        par_dist = math.sqrt(dist2) / sphere_rad

        if par_dist <= inner_par_rad:
            # The centre is atp_max.
            atp_val = max_atp
        else:
            # Discard the inner radius from our parametric distance.
            x0 = par_dist - inner_par_rad

            # Map distance to sigmoind domain min to max range.
            x1 = x0 * scaling_to_sigmoid_domain
            x2 = x1 + sigmoind_domain_min

            # Calculate sigmoid value and subtract the min sigmoid value for this domain.
            s = sigmoid(x2)
            s1 = s

            # Scale to atp range.
            atp_range = max_atp - atp_base
            atp_val = atp_base + s1 * atp_range

        newMapArray.SetValue(cid, atp_val)

sphere_rad = 1500
centre = (3000, 3022.5, 0)
max_atp = 0.25

# Do LO ATP 1 by running over all cell centress in polydata and processing only the ones that fall in the specified radius.
for cid in range (0, centres.GetNumberOfPoints()):

    # Distance from the point to sphere centre.
    dist2 = vtk.vtkMath.Distance2BetweenPoints(centres.GetPoint(cid), centre)

    # If the point falls into the sphere radius.
    if dist2 < sphere_rad**2:

        # Parametric distance from point to sphere centre.
        par_dist = math.sqrt(dist2) / sphere_rad

        if par_dist <= inner_par_rad:
            # The centre is atp_max.
            atp_val = max_atp
        else:
            # Discard the inner radius from our parametric distance.
            x0 = par_dist - inner_par_rad

            # Map distance to sigmoind domain min to max range.
            x1 = x0 * scaling_to_sigmoid_domain
            x2 = x1 + sigmoind_domain_min

            # Calculate sigmoid value and subtract the min sigmoid value for this domain.
            s = sigmoid(x2)
            s1 = s

            # Scale to atp range.
            atp_range = max_atp - atp_base
            atp_val = atp_base + s1 * atp_range

        newMapArray.SetValue(cid, atp_val)

sphere_rad = 1500
centre = (7000, 3022.5, 0)
max_atp = 0.3

# Do LO ATP 2 by running over all cell centres in polydata and processing only the ones that fall in the specified radius.
for cid in range (0, centres.GetNumberOfPoints()):

    # Distance from the point to sphere centre.
    dist2 = vtk.vtkMath.Distance2BetweenPoints(centres.GetPoint(cid), centre)

    # If the point falls into the sphere radius.
    if dist2 < sphere_rad**2:

        # Parametric distance from point to sphere centre.
        par_dist = math.sqrt(dist2) / sphere_rad

        if par_dist <= inner_par_rad:
            # The centre is atp_max.
            atp_val = max_atp
        else:
            # Discard the inner radius from our parametric distance.
            x0 = par_dist - inner_par_rad

            # Map distance to sigmoind domain min to max range.
            x1 = x0 * scaling_to_sigmoid_domain
            x2 = x1 + sigmoind_domain_min

            # Calculate sigmoid value and subtract the min sigmoid value for this domain.
            s = sigmoid(x2)
            s1 = s

            # Scale to atp range.
            atp_range = max_atp - atp_base
            atp_val = atp_base + s1 * atp_range

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

