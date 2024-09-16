#!/usr/bin/env python3

'''Script for computing sequence alignments using Needleman-Wunsch with
   linear gap penalties.
Arguments:
    f - FASTA file with sequences in FASTA format.
    s - JSON with the score matrix for alignment.
    d - The gap penalty for the alignment.

Outputs:
    Prints alignment to console.

Example Usage:
    python 2b.py -f sequences.fasta -s score_matrix.json -d 100
'''

import argparse
import json
import numpy as np


'''Computes the actual string alignments given the traceback matrix.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    t: the traceback matrix
Returns:
    a_x: the string for the alignment of x's sequence
    a_y: the string for the alignment of y's sequence
'''
def traceback(x, y, t):
    ''' Complete this function. '''
    DIAGONAL = 0
    LEFT = 1
    UP = 2
    a_x = ''
    a_y = ''
    i = len(x)
    j = len(y)
    while i > 0 or j > 0:
        if t[i][j] == DIAGONAL:
            a_x = x[i - 1] + a_x
            a_y = y[j - 1] + a_y
            i -= 1
            j -= 1
        elif t[i][j] == LEFT:
            a_x = x[i - 1] + a_x
            a_y = '-' + a_y
            i -= 1
        elif t[i][j] == UP:
            a_x = '-' + a_x
            a_y = y[j - 1] + a_y
            j -= 1
    return a_x, a_y
    


'''Computes the score and alignment of two strings.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    s: the score matrix
    d: the gap opening/extension penalty
Returns:
    score: the score of the optimal sequence alignment
    a_x: the aligned first string
    a_y: the aligned second string
The latter two are computed using the above traceback method.
'''
def sequence_alignment(x, y, s, d):
    # initializing
    DIAGONAL = 0
    LEFT = 1
    UP = 2
    
    ''' Recurrence matrix, redefine/use as necessary. '''
    m = np.zeros((len(x)+1, len(y)+1))
    ''' Traceback matrix, redefine/use as necessary. '''
    tb = np.zeros((len(x)+1, len(y)+1))

    for i in range(len(x)+1):
        m[i][0] = -d*i
        tb[i][0] = LEFT
    for j in range(len(y)+1):
        m[0][j] = -d*j
        tb[0][j] = UP

    # fills in from left to right, top to bottom
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            # index for the string starts at 0
            d_score = m[i-1][j-1] + s[x[i-1]][y[j-1]]
            l_score = m[i-1][j] - d
            u_score = m[i][j-1] - d
            m[i][j] = max(d_score, l_score, u_score)

            if max(d_score, l_score, u_score) == d_score:
                tb[i][j] = DIAGONAL
            elif max(d_score, l_score, u_score) == l_score:
                tb[i][j] = LEFT
            elif max(d_score, l_score, u_score) == u_score:
                tb[i][j] = UP

    a_x, a_y = traceback(x, y, tb)
    score = m[len(x)][len(y)]
    return score, (a_x, a_y)


'''Prints two aligned sequences formatted for convenient inspection.
Arguments:
    a_x: the first sequence aligned
    a_y: the second sequence aligned
Outputs:
    Prints aligned sequences (80 characters per line) to console
'''
def print_alignment(a_x, a_y):
    assert len(a_x) == len(a_y), "Sequence alignment lengths must be the same."
    for i in range(1 + (len(a_x) // 80)):
        start = i * 80
        end = (i + 1) * 80
        print(a_x[start:end])
        print(a_y[start:end])
        print()


def main():
    parser = argparse.ArgumentParser(
        description='Calculate sequence alignments for two sequences with a linear gap penalty.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-s', action="store", dest="s", type=str, required=True)
    parser.add_argument('-d', action="store", dest="d", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    score_matrix_file = args.s
    d = args.d

    with open(fasta_file) as f:
        _, x, _, y = [line.strip() for line in f.readlines()]
    with open(score_matrix_file) as f:
        s = json.loads(f.readlines()[0])

    score, (a_x, a_y) = sequence_alignment(x, y, s, d)
    print("Alignment:")
    print_alignment(a_x, a_y)
    print("Score: " + str(score))


if __name__ == "__main__":
    main()
