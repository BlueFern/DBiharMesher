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
    FlatMeshGenerator.xQuads = 16
    FlatMeshGenerator.yQuads = 16

    taskMesh = FlatMeshGenerator.buildMesh(1, 1, "quadMeshFullc256.vtp")

    ecMesh = FlatMeshGenerator.buildMesh(20, 4, "quadMeshFullECc256.vtp")
    smcMesh = FlatMeshGenerator.buildMesh(4, 52, "quadMeshFullSMCc256.vtp")

    atpMesh = FlatMeshGenerator.buildATPMesh(ecMesh, "quadMeshFullATPc256.vtp")

if __name__ == "__main__":
    main()

