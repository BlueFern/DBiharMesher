
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

/**
 * Show grids with optional centreline.
 */
void showGrids(std::vector<vtkSmartPointer<vtkStructuredGrid> > grids, vtkSmartPointer<vtkPolyData> centreline = NULL);

/**
 * Write polydata in XML format.
 */
void writePolyData(vtkPolyData *input, std::string fileName);
