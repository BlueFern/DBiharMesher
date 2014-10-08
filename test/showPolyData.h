class vtkPolyData;
class vtkStructuredGrid;

/**
 * Print an array to stdout as a 3D point.
 */
void PrintPoint(double *point);

/**
 * Show input and output polydata.
 */
void showPolyData(vtkPolyData *input, vtkStructuredGrid *output);
