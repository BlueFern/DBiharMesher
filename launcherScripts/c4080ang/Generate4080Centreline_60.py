# -*- coding: utf-8 -*-
"""
Generate centreline and write it out as .vtk legacy format.
"""

import os
import sys

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Import path for the CentrelineGenerator script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import CentrelineGenerator

# A centreline for a mesh with 4080 cores.
CentrelineGenerator.segmentList = [8.84,[(8.84,70),None,None],[(8.84,130),None,None]]
CentrelineGenerator.radiusBase = 1.2732395447351628
CentrelineGenerator.outputFileName = "c4080Centreline_60.vtk"
CentrelineGenerator.sphereRadius = None

def main():
    # CentrelineGenerator.GenerateCentreline(CentrelineGenerator.BuildDecreasingRadiiScalars)
    CentrelineGenerator.GenerateCentreline()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
