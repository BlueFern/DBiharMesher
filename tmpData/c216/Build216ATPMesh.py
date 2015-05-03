# -*- coding: utf-8 -*-
"""
Create an initial ATP Profile for an ECs Mesh and write it out as .vtp.
"""

import os
import sys

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Import path for the BuildATPMesh script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import BuildATPMesh

# This is for the c216 mesh.
meshFile = "quadMeshFullECc216.vtp"
atpFile = "quadMeshFullATPc216.vtp"
numBranches = 3
numQuads = 216
numECsPerCol = 4
atpGradient = 0.15
atpMin = 0.1
atpMax = 1.0

def main():
    BuildATPMesh.buildATPMesh()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
