# -*- coding: utf-8 -*-
"""
Write initial ATP Profile for an ECs Mesh in legacy VTK format as .vtk.
"""

import os
import sys
import shutil

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Relative import path for the DumpATPMeshToLegacyFormat script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import DumpATPMeshToLegacyFormat
reload(DumpATPMeshToLegacyFormat)

# This is for the c4080 mesh.
DumpATPMeshToLegacyFormat.numQuadsPerRing0 = 40

def main():
    for angle in range(60, 120, 10):
        print "Processing angle", angle
        os.chdir(str(angle))
        print "Current working directory:", os.getcwd()

        DumpATPMeshToLegacyFormat.taskMeshIn = 'quadMeshFullc4080.vtp'
        DumpATPMeshToLegacyFormat.ecMeshIn = 'quadMeshFullECc4080.vtp'
        DumpATPMeshToLegacyFormat.atpMeshIn = 'quadMeshFullATPc4080.vtp'

        DumpATPMeshToLegacyFormat.writeATPLegacyVTK()

        os.chdir('..')

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
else:
    print __file__, "is to be run as main script."
