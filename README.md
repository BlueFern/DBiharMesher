DbiharPatchFilter
=================

vtkDbiharPatchFilter can be used for generating quadrilateral meshes for surface patches. Surface patches can be combined to form surface meshes for complex geometric objects.

vtkDbiharPatchFilter is a [VTK](http://vtk.org) filter/wrapper for biharmonic equation solver library Dbihar written in Fortran. Dbihar Fortran code was downloaded from [Netlib](http://www.netlib.org/bihar/index.html). The Dbihar library solves the biharmonic equations in the parametric space.

The filer is intended to simplify the use of Dbihar library in the in the process of generating quadrilateral meshes.

Input
-----

The input to vtkDbiharPatchFilter is a vtkPolyData object which contains a vtkPolyLine over a set of points contained within the vtkPolyData object. In addition to the vtkPolyLine boundary the shape of the output patch also depends on the spatial derivatives associated with the vtkPolyData points.

Output
------

The output from the vtkDbiharPatchFilter filter is a set of nodes for a structured grid representing the patch surface. The set of nodes is converted to vtkStructuredGrid for visualisation.

