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

# This is for the c4000 mesh.
GenerateATPMapV2.centrelineFile = "c512Centreline.vtk"
GenerateATPMapV2.meshFile = "quadMeshFullECc512.vtp"
GenerateATPMapV2.atpFile = "quadMeshFullATPc512.vtp"
GenerateATPMapV2.debugAtpFile = "quadMeshFullATPV2Debugc512.vtp"
GenerateATPMapV2.numBranches = 1
GenerateATPMapV2.numQuads = 512
GenerateATPMapV2.numAxialQuads = 26
GenerateATPMapV2.numECsPerCol = 4
GenerateATPMapV2.atpGradient = 3.0
GenerateATPMapV2.atpMin = 0.3
GenerateATPMapV2.atpMax = 1.0

def main():
    GenerateATPMapV2.buildATPMesh()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
