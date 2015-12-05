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

    FlatMeshGenerator.xQuads = 40
    FlatMeshGenerator.yQuads = 34
        
    taskMesh = FlatMeshGenerator.buildMesh(1, 1, "quadMeshFullc4080.vtp")
    
    ecMesh = FlatMeshGenerator.buildMesh(20, 4, "quadMeshFullECc4080.vtp")
    smcMesh = FlatMeshGenerator.buildMesh(4, 52, "quadMeshFullSMCc4080.vtp")

    atpMesh = FlatMeshGenerator.buildATPMesh(ecMesh, "quadMeshFullATPc4080.vtp")
    
if __name__ == "__main__":
    main()
    
