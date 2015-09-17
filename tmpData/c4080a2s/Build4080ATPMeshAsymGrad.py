# -*- coding: utf-8 -*-
"""
Create an initial ATP Profile for an ECs Mesh and write it out as .vtp.
"""

import os
import sys

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Import path for the BuildATPMeshA2 script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import BuildATPMeshA2

# This is for the c4080 mesh.
BuildATPMeshA2.meshFile = "quadMeshFullECc4080.vtp"
BuildATPMeshA2.atpFile = "quadMeshFullATPc4080.vtp"
BuildATPMeshA2.numBranches = 3
BuildATPMeshA2.numQuads = 4080
BuildATPMeshA2.numECsPerCol = 4
BuildATPMeshA2.atpGradient0 = 0.02
BuildATPMeshA2.atpGradient1 = 0.02
BuildATPMeshA2.atpGradient2 = 0.006666666666666666
BuildATPMeshA2.atpMin = 0.1
BuildATPMeshA2.atpMax = 1.0

def main():
    BuildATPMeshA2.buildATPMeshAsymGrad()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
