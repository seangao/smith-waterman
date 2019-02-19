#!/usr/bin/python
__author__ = "Sean Gao"
__email__ = "sean.gao@yale.edu"
__copyright__ = "Copyright 2019"
__license__ = "GPL"
__version__ = "1.0.0"
### Usage:    smith_waterman -i <input file> -s <score file>
### Example:  smith_waterman -i input.txt -s blosum62.txt
### Note:     Smith-Waterman Algorithm

import argparse
import numpy as np

### Implement Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap=-2, extGap=-1):
    ### Parse sequences
    seqs = []
    with open(inputFile) as f:
        seqs = f.readlines()
    seqs = seqs + ['']
    seqs = [line.strip().upper() for line in seqs if line.strip()]

    ### Convert score file into matrix
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

    ### Align matrix setup
    a = seqs[0]
    b = seqs[1]
    n = len(a)
    m = len(b)
    H = np.zeros((n + 1, m + 1))
    state = np.zeros((n + 1, m + 1))

    ### Create gap penalty vector
    W = np.arange(openGap, openGap + extGap * max(n, m), extGap)

    ### Create align matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = H[i - 1][j - 1] + S[SRows[a[i - 1]], SCols[b[j - 1]]]
            insertion = max([H[i - k][j] + W[k - 1] for k in range(1, i + 1)])
            deletion = max([H[i][j - l] + W[l - 1] for l in range(1, j + 1)])
            H[i][j] = max(match,
                insertion,
                deletion,
                0)
            if H[i][j] == match:
                state[i][j] = 1
            elif H[i][j] == insertion:
                state[i][j] = 2
            elif H[i][j] == deletion:
                state[i][j] = 3

    ### Traceback
    ind = np.unravel_index(np.argmax(H, axis=None), H.shape)
    score = H[ind]
    anew = a[ind[0] - 1]
    bnew = b[ind[1] - 1]
    while H[ind] != 0 and ind[0] > 1 and ind[1] > 1:
        if state[ind] == 1:
            ind = ind[0] - 1, ind[1] - 1
            anew += a[ind[0] - 1]
            bnew += b[ind[1] - 1]
        elif state[ind] == 2:
            ind = ind[0] - 1, ind[1]
            anew += a[ind[0] - 1]
            bnew = bnew[:-1] + '-' + bnew[-1]
        elif state[ind] == 3:
            ind = ind[0], ind[1] - 1
            anew = anew[:-1] + '-' + anew[-1]
            bnew += b[ind[1] - 1]
        else:
            break

    bnew = bnew[::-1]
    anew = anew[::-1]
    align = ''
    for i in range(len(anew)):
        if anew[i] == bnew[i]:
            align += '|'
        else:
            align += ' '

    print (" ---------------------------------------------------")
    print ("|       Sequences                                   |")
    print (" ---------------------------------------------------")
    print ("")
    print ("sequence1\t%s" % a)
    print ("sequence2\t%s" % b)
    print ("")

    print (" ---------------------------------------------------")
    print ("|       Score Matrix                                |")
    print (" ---------------------------------------------------")
    print ("")
    print ('\t'.join(['', ''] + list(b)))
    a = [''] + list(a)
    for i in range(H.shape[0]):
        print ('\t'.join([a[i]] + [str(k) for k in H[i].tolist()]))
    print ("")

    print (" ---------------------------------------------------")
    print ("|       Best Local Alignment                        |")
    print (" ---------------------------------------------------")
    print ("")
    print ("Alignment Score:\t\t%s" % score)
    print (anew)
    print (align)
    print (bnew)

def main():
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

    ### Print input and score file names.
    print ("input file : %s" % inputFile)
    print ("score file : %s" % scoreFile)
    print ("open gap penalty : %s" % openGap)
    print ("extension gap penalty : %s" % extGap)

    ### Run Smith-Waterman Algorithm
    runSW(args.input, args.score, args.opengap, args.extgap)

if __name__ == '__main__':
    main()
