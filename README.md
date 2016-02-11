# DbiharPatchFilter


vtkDbiharPatchFilter can be used for generating quadrilateral meshes for surface patches. Surface patches can be combined to form surface meshes for complex geometric objects.

vtkDbiharPatchFilter is a [VTK](http://vtk.org) filter/wrapper for biharmonic equation solver library Dbihar written in Fortran. Dbihar Fortran code was downloaded from [Netlib](http://www.netlib.org/bihar/index.html). The Dbihar library solves the biharmonic equations in the parametric space.

The filer is intended to simplify the use of Dbihar library in the in the process of generating quadrilateral meshes.

## Input

The input to vtkDbiharPatchFilter is a vtkPolyData object which contains a vtkPolyLine over a set of points contained within the vtkPolyData object. In addition to the vtkPolyLine boundary the shape of the output patch also depends on the spatial derivatives associated with the vtkPolyData points.

## Output

The output from the vtkDbiharPatchFilter filter is a set of nodes for a structured grid representing the patch surface. The set of nodes is converted to vtkStructuredGrid for visualisation.


#DbiharPatchSmooth


vtkDbiharPatchSmooth is a filter, generated to smooth the arrangement of quads of an simplyfied coronary artery bifurcation model (such as quadMeshFullc4080.vtp). To do so the algorithm must be supplied with the original model as input. It then generates two new daugther branches, appends them to the existing trunk and gives that new model as output. It creates new daughter branches by taking their contours and applying the DbiharPatchFilter on each of them. 
To be able to calculate the necessary data from the original model, the number of radial quads ('numRadialQuads', halved) must be set to the same value as it is in the pipeline that creates the input data (i.e. meshGen4080.cxx). With the aim of this parameter the algorithm is able to calculate all the necessary data by itself. First it gets important coordinates (point in the centre of bifurcation and endpoints of daughter branches) and calculates the models properties (lengths, radii, angles). Knowing that there are four endothelial cells per quad in axial direction and the lengths of the branches, the number of axial quads per branch ('numAxialQuads0', 'numAxialQuads1', 'numAxialQuads2') can be determined. In order to define the contour of a single daugther branch, consecutive points are appended into an array. This happens by applying specific functions on the original point array (the original model



![Image name/tag goes here](https://raw.github.com/bluefern/dbiharmesher/master/doc/images/boundaryWithDerivatives_numbered.png)
