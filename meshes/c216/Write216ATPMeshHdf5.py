# -*- coding: utf-8 -*-
"""
Write initial ATP Profile for an ECs Mesh in legacy VTK format as .vtk.
"""

import os
import sys

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Relative import path for the DumpATPMeshToLegacyFormat script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import DumpATPToHdf5

# This is for the flat c24 mesh.

DumpATPToHdf5.axialQuads = 6
DumpATPToHdf5.circQuads = 12

DumpATPToHdf5.numECsPerCol = 4
DumpATPToHdf5.numSMCsPerRow = 4
DumpATPToHdf5.taskMeshIn = "quadMeshFullc216.vtp"
DumpATPToHdf5.ecMeshIn = "quadMeshFullECc216.vtp"
DumpATPToHdf5.atpMeshIn = "quadMeshFullATPc216.vtp"

def main():
    DumpATPToHdf5.writeHdf5()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
else:
    print __file__, "is to be run as main script."
