# -*- coding: utf-8 -*-
"""
Create an initial ATP Profile for an ECs Mesh and write it out as .vtp.
"""

import os
import sys

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Import path for the GenerateATPMapV2 script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import GenerateATPMapV2

# This is for the c216 mesh.
GenerateATPMapV2.centrelineFile = "c216Centreline.vtk"
GenerateATPMapV2.meshFile = "quadMeshFullECc216.vtp"
GenerateATPMapV2.atpFile = "quadMeshFullATPc216.vtp"
GenerateATPMapV2.debugAtpFile = "quadMeshFullATPV2Debugc216.vtp"
GenerateATPMapV2.numBranches = 3
GenerateATPMapV2.numQuads = 216
GenerateATPMapV2.numAxialQuads = 6
GenerateATPMapV2.numECsPerCol = 4
GenerateATPMapV2.atpGradient = 4.0
GenerateATPMapV2.atpMin = 0.1
GenerateATPMapV2.atpMax = 1.0

def main():
    GenerateATPMapV2.buildATPMesh()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
