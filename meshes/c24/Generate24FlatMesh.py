import os
import sys

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Import path for the CentrelineGenerator script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import FlatMeshGenerator

def main():
        
    FlatMeshGenerator.buildMesh(1, 1, "quadMeshFullc24.vtp")
    polydata = FlatMeshGenerator.buildMesh(20, 4, "quadMeshFullECc24.vtp")
    FlatMeshGenerator.buildMesh(4, 52, "quadMeshFullSMCc24.vtp")
    FlatMeshGenerator.buildATPMesh(polydata, "quadMeshFullATPc24.vtp")
    

if __name__ == "__main__":
    main()
    
