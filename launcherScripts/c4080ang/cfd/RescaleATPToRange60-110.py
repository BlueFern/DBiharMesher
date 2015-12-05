# -*- coding: utf-8 -*-
"""
Create initial ATP profile based on a CDF simulation.
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

import RescaleATPToRange
reload(RescaleATPToRange)

RescaleATPToRange.outMin = 0.2
RescaleATPToRange.outMax = 0.8

def main():
    for angle in range(60, 120, 10):
        RescaleATPToRange.inputFile = os.path.join(str(angle), 'quadMeshFullATPc4080_' + str(angle) + '.vtp')
        print 'Input file:', RescaleATPToRange.inputFile

        RescaleATPToRange.inputSurfaceFile = os.path.join(str(angle), 'quadMeshFullECc4080.vtp')
        print 'Input surface file:', RescaleATPToRange.inputSurfaceFile

        RescaleATPToRange.outputFile = os.path.join(str(angle), 'quadMeshFullATPc4080.vtp')
        print 'Output file:', RescaleATPToRange.outputFile

        RescaleATPToRange.outputSurfaceFile = 'quadMeshFullATPc4080_' + str(angle) + '.vtp'
        print 'Output surface file:', RescaleATPToRange.outputSurfaceFile

        RescaleATPToRange.rescaleATPToRange()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
else:
    print __file__, "is to be run as main script."