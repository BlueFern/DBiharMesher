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

# Large centreline.
CentrelineGenerator.radiusBase = 2.5
CentrelineGenerator.segmentList = [20,
               [(20,60),[(20,15),[(30,350),None,None],[(30,60),None,None]],[(40,105),None,None]],
               [(20,150),[(50,105),None,None],[(30,150),[(40,105),None,None],[(60,150),None,None]]]]
               
CentrelineGenerator.outputFileName = "decreasingRadiiCentreline.vtk"
CentrelineGenerator.sphereRadius = None

def main():
    CentrelineGenerator.GenerateCentreline(CentrelineGenerator.BuildDecreasingRadiiScalars)

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)