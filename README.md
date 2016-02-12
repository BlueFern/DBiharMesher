# DbiharPatchFilter


vtkDbiharPatchFilter can be used for generating quadrilateral meshes for surface patches. Surface patches can be combined to form surface meshes for complex geometric objects.

vtkDbiharPatchFilter is a [VTK](http://vtk.org) filter/wrapper for biharmonic equation solver library Dbihar written in Fortran. Dbihar Fortran code was downloaded from [Netlib](http://www.netlib.org/bihar/index.html). The Dbihar library solves the biharmonic equations in the parametric space.

The filer is intended to simplify the use of Dbihar library in the in the process of generating quadrilateral meshes.

## Input

The input to vtkDbiharPatchFilter is a vtkPolyData object which contains a vtkPolyLine over a set of points contained within the vtkPolyData object. In addition to the vtkPolyLine boundary the shape of the output patch also depends on the spatial derivatives associated with the vtkPolyData points.

## Output

The output from the vtkDbiharPatchFilter is a set of nodes for a structured grid representing the patch surface. The set of nodes is converted to vtkStructuredGrid for visualisation.


# Generation of a Bifurcation Model

*The examples in this chapter are from a model with 4080 quads*


## Centreline

The centreline for the bifurcation model is created by the Python program 'Generate4080Centreline.py'. It contains the branches lengths, angles and radius as input data and passes them to 'CentrelineGenerator.py'. **TODO: difference consistent/decreasing radii (Murray's Law)** The output centreline is two dimensional and the radii in every point are set as scalars. The resulting vtkPolyData is then saved in a file.
 
![Centreline](https://raw.github.com/bluefern/dbiharmesher/master/doc/images/centreline_4080.png "Centreline")
 
## Mesh Generation

The mesh generator program 'meshGen4080.cxx' starts by reading in the centreline file ('c4080Centreline.vtk') and setting the number of radial quads. It prepares the data for the filter that actually generates the surface points, the 'dbiharPatchFilter' (*vtkCentrelineToDbiharPatch*?). The points are then used to form the quads.
This pipeline is a series of scripts and its input can be altered to get different sized and shaped bifurcation models. 

![Old_Mesh](https://raw.github.com/bluefern/dbiharmesher/master/doc/images/old_Mesh.png "Example for a first solution mesh")


## Mesh Smoothing

The solution from the mesh generator can be improved with the optional filter 'vtkDbiharPatchSmooth'. It has been developed after undesired distortions appeared in simulations run on the bifurcation model. Some of the joined quads form sharp angles along the daughter branches towards the centre of the bifurcation and the aim of this filter is to smoothen them. To do so it reads in the vtkPolyData solution (i.e. 'quadMeshFullc4080.vtp'), generates two new daughter branches on the basis of this input, attaches them to the existing parent branch and gives the new model as output. The lines of the mesh are more consitent due to the way the boundaries for each branch is defined before the 'vtkDbiharPatchFilter' is applied. The boundary line itself is defined as a consecutive collection of points from the temporary solution (output of meshGenerator) and surrounds one single branch (see following image "Boundary"). Each point is then equipped with a derivative. The magnitude and direction of these vectors are functions and dependent on the properties of the input model. They affect the propagations of the generated lines as well as the shape of the solution. 

![consecutiveOrderOfPoints](https://raw.github.com/bluefern/dbiharmesher/master/doc/images/boundaryWithDerivatives_numbered.png "Boundary")

![New_Mesh](https://raw.github.com/bluefern/dbiharmesher/master/doc/images/new_Mesh.png "Example for improved mesh")

**TODO: Other Meshes (SMCs, ECs) and (rectangular, brick lay)!?**
