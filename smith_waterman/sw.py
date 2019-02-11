#!/usr/bin/python
__author__ = "Sean Gao"
__email__ = "sean.gao@yale.edu"
__copyright__ = "Copyright 2019"
__license__ = "GPL"
__version__ = "1.0.0"
### Usage:    python sw.py -i <input file> -s <score file>
### Example:  python sw.py -i input.txt -s blosum62.txt
### Note:     Smith-Waterman Algorithm

import argparse
import numpy as np

### Read input file and score file.
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()

### Implement Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
### Print input and score file names.
print ("input file : %s" % inputFile)
print ("score file : %s" % scoreFile)
print ("open gap penalty : %s" % openGap)
print ("extension gap penalty : %s" % extGap)
### Run Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)
