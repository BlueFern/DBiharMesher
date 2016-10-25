# -*- coding: utf-8 -*-
"""
Write multiple ATP profiles to hdf5 format.
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

import MultipleATPToHDF5

# This is for the flat c4080 mesh.

MultipleATPToHDF5.axialQuads = 34
MultipleATPToHDF5.circQuads = 40
MultipleATPToHDF5.timeStep = 0.1
MultipleATPToHDF5.numECsPerCol = 4
MultipleATPToHDF5.numSMCsPerRow = 4
MultipleATPToHDF5.taskMeshIn = "quadMeshFullc4080.vtp"



def main():

    MultipleATPToHDF5.atpMeshPattern = "50_new/quadMeshFullECc4080_*.vtp"
    MultipleATPToHDF5.writeHdf5()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
else:
    print __file__, "is to be run as main script."
