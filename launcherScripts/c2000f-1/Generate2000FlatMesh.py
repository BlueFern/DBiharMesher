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

    FlatMeshGenerator.bifurcation = False
    FlatMeshGenerator.xQuads = 50
    FlatMeshGenerator.yQuads = 40
        
    taskMesh = FlatMeshGenerator.buildMesh(1, 1, "quadMeshFullc2000.vtp")
    
    ecMesh = FlatMeshGenerator.buildMesh(20, 4, "quadMeshFullECc2000.vtp")
    smcMesh = FlatMeshGenerator.buildMesh(4, 52, "quadMeshFullSMCc2000.vtp")

    atpMesh = FlatMeshGenerator.buildATPMesh(ecMesh, "quadMeshFullATPc2000.vtp")
    
if __name__ == "__main__":
    main()
    
