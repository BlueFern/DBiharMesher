# -*- coding: utf-8 -*-
"""
Generate centreline and write it out as .vtk legacy format.
"""

import os
import sys

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Import path for the CentrelineGenerator script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import CentrelineGenerator

# A centreline for a mesh with 216 cores.
CentrelineGenerator.segmentList = [1.7,[(1.7,60),None,None],[(1.7,120),None,None]]
CentrelineGenerator.radiusBase = 0.382
CentrelineGenerator.outputFileName = "c216Centreline.vtk"
CentrelineGenerator.sphereRadius = None

def main():
    # CentrelineGenerator.GenerateCentreline(CentrelineGenerator.BuildDecreasingRadiiScalars)
    CentrelineGenerator.GenerateCentreline()

if __name__ == '__main__':
    print("Starting", os.path.basename(__file__))
    main()
    print("Exiting", os.path.basename(__file__))
