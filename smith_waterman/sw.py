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

inputFile = args.input
scoreFile = args.score
openGap = args.opengap
extGap = args.extgap

### Implement Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
    seqs = []
    with open(inputFile) as f:
        seqs = f.readlines()
    seqs = seqs + ['']
    seqs = [line.strip() for line in seqs if line.strip()]

    scoreRows = []
    scoreCols = []
    S = []
    with open(scoreFile) as f:
        lines = f.readlines()
        scoreCols = lines[0].strip().split()
        for line in lines[1:]:
            line = line.strip().split()
            if line:
                scoreRows.append(line[0])
                row = [float(i) for i in line[1:]]
                S.append(row)
    SRows = {k: v for v, k in enumerate(scoreRows)}
    SCols = {k: v for v, k in enumerate(scoreCols)}
    S = np.matrix(S)

    a = seqs[0]
    b = seqs[1]
    n = len(a)
    m = len(b)
    H = np.zeros((n + 1, m + 1))

    W = np.arange(openGap, openGap + extGap * max(n, m), extGap)

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            H[i][j] = max(H[i - 1][j - 1] + S[SRows[a[i - 1]], SCols[b[j - 1]]],
                max([H[i - k][j] + W[k - 1] for k in range(1, i + 1)]),
                max([H[i][j - l] + W[l - 1] for l in range(1, j + 1)]),
                0)

    ind = np.unravel_index(np.argmax(H, axis=None), H.shape)
    np.set_printoptions(threshold=np.nan)
    print H
    print H[ind]
    anew = a[ind[0] - 1]
    bnew = b[ind[1] - 1]
    while H[ind] != 0:
        next = max(H[ind[0] - 1, ind[1] - 1],
            H[ind[0] - 1, ind[1]],
            H[ind[0], ind[1] - 1])
        if next == H[ind[0] - 1, ind[1] - 1]:
            ind = ind[0] - 1, ind[1] - 1
            anew += a[ind[0] - 1]
            bnew += b[ind[1] - 1]
        elif next == H[ind[0] - 1, ind[1]]:
            ind = ind[0] - 1, ind[1]
            anew += a[ind[0] - 1]
            bnew = bnew[:-1] + '-' + bnew[-1]
        else:
            ind = ind[0], ind[1] - 1
            anew = anew[:-1] + '-' + anew[-1]
            bnew += b[ind[1] - 1]

    print anew[::-1]
    print bnew[::-1]



### Print input and score file names.
print ("input file : %s" % inputFile)
print ("score file : %s" % scoreFile)
print ("open gap penalty : %s" % openGap)
print ("extension gap penalty : %s" % extGap)
### Run Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)
