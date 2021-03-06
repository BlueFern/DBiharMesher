# -*- coding: utf-8 -*-
"""
Write EC and SMC Meshes in legacy VTK format as .vtk.
"""

import os
import sys

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Relative import path for the DumpMeshToLegacyFormat script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import DumpMeshToLegacyFormat

# This is for a flat c24 mesh.
DumpMeshToLegacyFormat.numQuadsPerRing = 4
DumpMeshToLegacyFormat.meshSet = [
"quadMeshFullc24.vtp",
"quadMeshFullECc24.vtp",
"quadMeshFullSMCc24.vtp"
]

def main():
    DumpMeshToLegacyFormat.writeLegacyVTK()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
else:
    print __file__, "is to be run as main script."
