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

RescaleATPToRange.inputFile = 'quadMeshFullJPLCc4080FlowRate050.vtp'
RescaleATPToRange.outputFile = 'quadMeshFullATPc4080.vtp'
RescaleATPToRange.outMin = 0.2
RescaleATPToRange.outMax = 1.0

def main():
    RescaleATPToRange.rescaleATPToRange()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)