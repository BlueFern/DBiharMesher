
class vtkPolyData;
class vtkStructuredGrid;

/**
 * Print an array to stdout as a 3D point.
 */
void PrintPoint(double *point);

/**
 * Show input and output polydata.
 */
void showPolyData(vtkPolyData *input, vtkStructuredGrid *output, double derivateScaling = 0.1);

/**
 * Show polydata with vectors.
 */
void showPolyData1(vtkPolyData *input, double vectorScaling = 1.0);
