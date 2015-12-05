# -*- coding: utf-8 -*-
"""
For each saddle point (out of three) in a bifurcation EC mesh place an 'imaginary' sphere around the point.
Iterate over all EC centres and if a cell centre falls into one of the spheres, map the distance to the
sphere centre to an ATP value.
"""

import os
import vtk
import math

import numpy as np
import matplotlib.pyplot as plt

numECsCirc = 20
numECsAx = 4

atp_inc = 0.02
atp_base = 0.35
atp_max = [0.25, 0.3, 0.68] # 0, 1, 2.

sigmoid_grad = -3.0
sphere_rads = [1900, 1900, 1000]
rad_inc = 0.15
rad_dec = 0.06
inner_rad = 1/2.0
sigmoind_domain_min = -2.0
sigmoind_domain_max = 2.0
sigmoid_domain = sigmoind_domain_max - sigmoind_domain_min
scaling_to_sigmoid_domain = sigmoid_domain / (1.0 - inner_rad)

# This is for the c4080 mesh.
numQuadRings = 34
numQuadsPerRing = 40

def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-sigmoid_grad * x))

def main():

    angle_num = 1

    for angle in range(60, 120, 10):

        print 'Processing angle', angle, '...'
        # Read the EC mesh.
        ecFileName = os.path.join(str(angle), 'quadMeshFullECc4080.vtp')
        ecMeshReader = vtk.vtkXMLPolyDataReader()
        ecMeshReader.SetFileName(ecFileName)
        ecMeshReader.Update()

        # Remember the EC mesh.
        ecMesh = ecMeshReader.GetOutput()

        # Get EC centres.
        cellCentresFilter = vtk.vtkCellCenters()
        cellCentresFilter.VertexCellsOn()
        cellCentresFilter.SetInput(ecMesh)
        cellCentresFilter.Update()

        ecMeshCentres = cellCentresFilter.GetOutput()
        pointsPerBranch = ecMeshCentres.GetNumberOfPoints() / 3

        # Find three saddle points.
        # The first point is 3/4th of the EC row away from the end of the parent branch.
        # The second point is 1/4th of the EC row away from the end of the parent branch.
        # The third point is 3/4th of the way from the start of the first sibling branch.
        saddle_ids = [pointsPerBranch - (numQuadsPerRing / 4) * (numECsAx * numECsCirc),
                         pointsPerBranch - ((numQuadsPerRing / 4) * 3) * (numECsAx * numECsCirc),
                         pointsPerBranch + ((numQuadsPerRing / 4) * 3) * (numECsAx * numECsCirc)]

        saddle_points = [ecMeshCentres.GetPoint(sId) for sId in saddle_ids]

        atpArray = vtk.vtkDoubleArray()
        atpArray.SetName("initialATP")

        # Loop through each centre point.
        for cId in range(pointsPerBranch * 3):
            # If it falls within one of the three spheres, map it to ATP.
            point = ecMeshCentres.GetPoints().GetPoint(cId)

            # Find closest saddle point to work with.
            # List of distances to all three saddle points.
            dist2_list = [vtk.vtkMath.Distance2BetweenPoints(point, saddle_point) for saddle_point in saddle_points]

            # Position of the closest saddle point.
            saddle_id = dist2_list.index(min(dist2_list))

            sphere_rad = sphere_rads[saddle_id]
            max_atp = atp_max[saddle_id]

            # Decrement radius if on the first one; this is the changing angle for bifurcation outer curve.
            if saddle_id == 0:
                sphere_rad -= (angle_num * rad_dec) * sphere_rad

            # Increment radius if on the third one; this is the neck of the bifurcation.
            if saddle_id == 2:
                sphere_rad += (angle_num * rad_inc) * sphere_rad

            # Increase intensity if on the third one.
            if saddle_id == 2:
                max_atp += ((angle_num - 1) * atp_inc)

            # Does dist2 fall into sphere_rad2?
            if dist2_list[saddle_id] <= sphere_rad**2:

                # Map to parametric distance in the range [0, 1].
                par_dist = math.sqrt(dist2_list[saddle_id]) / sphere_rad

                # Map to ATP range, depending on the saddle id.
                if par_dist <= inner_rad:
                    # The centre is atp_max.
                    atp_val = max_atp
                else:
                    # Discard the inner radius from our parametric distance.
                    x0 = par_dist - inner_rad

                    # Map distance to sigmoind domain min to max range.
                    x1 = x0 * scaling_to_sigmoid_domain
                    x2 = x1 + sigmoind_domain_min

                    # Calculate sigmoid value and subtract the min sigmoid value for this domain.
                    s = sigmoid(x2)
                    s1 = s

                    # Scale to atp range.
                    atp_range = max_atp - atp_base
                    atp_val = atp_base + s1 * atp_range

            # Outside the sphere_rad2.
            else:
                atp_val = atp_base

            # Remember the ATP value.
            atpArray.InsertNextValue(atp_val)

        ecMeshCentres.GetCellData().SetScalars(atpArray)
        angle_num += 1

        ecMesh.GetCellData().SetScalars(atpArray)

        atpFileName = os.path.join(str(angle), 'quadMeshFullATPc4080.vtp')
        # atpFileName = 'quadMeshFullATPc4080_' + str(angle) + '.vtp'
        print 'Writing', atpFileName, '...'

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(atpFileName)
        writer.SetInput(ecMeshCentres)
        # writer.SetInput(ecMesh)
        writer.Update()

        # # Setup actor and mapper
        # mapper = vtk.vtkPolyDataMapper()
        # if vtk.VTK_MAJOR_VERSION <= 5:
        #     mapper.SetInput(ecMesh)
        # else:
        #     mapper.SetInputData(ecMesh)
        #
        # actor = vtk.vtkActor()
        # actor.SetMapper(mapper)
        #
        # # Setup render window, renderer, and interactor
        # renderer = vtk.vtkRenderer()
        # renderWindow = vtk.vtkRenderWindow()
        # renderWindow.AddRenderer(renderer)
        # renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        # renderWindowInteractor.SetRenderWindow(renderWindow)
        # renderer.AddActor(actor)
        # renderWindow.Render()
        # renderWindowInteractor.Start()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)

    #t1 = np.arange(sigmoind_domain_min, sigmoind_domain_max, 0.001)
    #plt.plot(t1, sigmoid(t1))
    #plt.show()

    main()

    print "Exiting", os.path.basename(__file__)
else:
    print __file__, "is to be run as main script."
