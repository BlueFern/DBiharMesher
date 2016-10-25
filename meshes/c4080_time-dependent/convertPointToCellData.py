# -*- coding: utf-8 -*-
"""
Create initial ATP profile based on a CDF simulation.
"""

import os
import sys

# Run in current directory.
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Import path for the BuildATPMesh script.
importPath = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../util'))
if not importPath in sys.path:
    sys.path.insert(1, importPath)
del importPath

import pointToCellData
reload(pointToCellData)


def main():

	pointToCellData.angle = 50
        pointToCellData.pointToCellData()

if __name__ == '__main__':
    print "Starting", os.path.basename(__file__)
    main()
    print "Exiting", os.path.basename(__file__)
else:
    print __file__, "is to be run as main script."
