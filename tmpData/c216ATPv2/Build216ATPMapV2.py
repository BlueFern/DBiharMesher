# -*- coding: utf-8 -*-
"""
Create an initial ATP Profile for an ECs Mesh and write it out as .vtp.
"""

import os
import sys

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Import path for the GenerateATPMap script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import GenerateATPMap

# This is for the c216 mesh.
GenerateATPMap.centrelineFile = "c216Centreline.vtk"
GenerateATPMap.meshFile = "quadMeshFullECc216.vtp"
GenerateATPMap.debugAtpFile = "quadMeshFullATPV2c216.vtp"
GenerateATPMap.atpFile = "quadMeshFullATPc216.vtp"
GenerateATPMap.numBranches = 3
GenerateATPMap.numQuads = 216
GenerateATPMap.numAxialQuads = 6
GenerateATPMap.numECsPerCol = 4
GenerateATPMap.atpGradient = 0.15
GenerateATPMap.atpMin = 0.1
GenerateATPMap.atpMax = 1.0

def main():
    GenerateATPMap.buildATPMesh()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
