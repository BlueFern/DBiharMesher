#include <stdlib.h>
#include <iostream>

int main(int argc, char *argv[])
{

	// TODO: This is to be extended to work with submeshes (ECs and SMCs meshes).

	if(argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " " << "<inputMeshFileName>" << std::endl;
		return EXIT_FAILURE;
	}

	// Read the mesh.

	// Reorder the quads.

	// Save these files in ASCII txt format. Sigh...
#if 0
"files/parent_points.txt");
"files/parent_cells.txt");

"files/left_daughter_points.txt");
"files/left_daughter_cells.txt");

"files/right_daughter_points.txt");
"files/right_daughter_cells.txt");


"files/parent_smc_mesh_points.txt");
"files/parent_smc_mesh_cells.txt");

"files/left_daughter_smc_mesh_points.txt");
"files/left_daughter_smc_mesh_cells.txt");

"files/right_daughter_smc_mesh_points.txt");
"files/right_daughter_smc_mesh_cells.txt");


"files/parent_ec_mesh_points.txt");
"files/parent_ec_mesh_cells.txt");

"files/left_daughter_ec_mesh_points.txt");
"files/left_daughter_ec_mesh_cells.txt");

"files/right_daughter_ec_mesh_points.txt");
"files/right_daughter_ec_mesh_cells.txt");


"files/parent_ec_centeroid_points.txt");
"files/parent_ec_centeroid_cells.txt");

"files/left_daughter_ec_centeroid_points.txt");
"files/left_daughter_ec_centeroid_cells.txt");

"files/right_daughter_ec_centeroid_points.txt");
"files/right_daughter_ec_centeroid_cells.txt");
#endif

	return EXIT_SUCCESS;
}
