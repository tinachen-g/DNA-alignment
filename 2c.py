#!/usr/bin/env python3

'''Script for computing sequence alignments using Needleman-Wunsch with
   affine gap penalties.
Arguments:
    f - FASTA file with sequences in FASTA format.
    s - JSON with the score matrix for alignment.
    d - The gap opening penalty for the alignment.
    e - The gap extension penalty for the alignment.

Outputs:
    Prints alignment to console.

Example Usage:
    python 2c.py -f sequences.fasta -s score_matrix.json -d 430 -e 30
'''

import argparse
import json
import numpy as np


'''Computes the actual string alignments given the traceback matrix.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    t: the traceback matrix, which stores values that point to which
       prior matrix was used to reach a given location in each of the
       3 matrices.
Returns:
    a_x: the string for the alignment of x's sequence
    a_y: the string for the alignment of y's sequence
'''
def traceback(x, y, t):
    ''' Complete this function. '''
    MID = 0
    BOT = 1
    TOP = 2
    a_x = ''
    a_y = ''
    i = len(x)
    j = len(y)
    curr_layer = int(t[i][j][MID])
    
    while i > 0 or j > 0:
        if curr_layer == MID:
            curr_layer = int(t[i][j][curr_layer])
            a_x = x[i - 1] + a_x
            a_y = y[j - 1] + a_y
            i -= 1
            j -= 1
            
        elif curr_layer == BOT:
            curr_layer = int(t[i][j][curr_layer])
            a_x = x[i - 1] + a_x
            a_y = '-' + a_y
            i -= 1

        elif curr_layer == TOP:
            curr_layer = int(t[i][j][curr_layer])
            a_x = '-' + a_x
            a_y = y[j - 1] + a_y
            j -= 1
        else: 
            raise Exception("Not possible path")

    return a_x, a_y


'''Computes the score and alignment of two strings using an affine gap penalty.
Arguments:
    x: the first string we're aligning
    y: the second string we're aligning
    s: the score matrix
    d: the gap opening penalty
    e: the gap extension penalty
Returns:
    score: the score of the optimal sequence alignment
    a_x: the aligned first string
    a_y: the aligned second string
The latter two are computed using the above traceback method.
'''
def affine_sequence_alignment(x, y, s, d, e):
    ''' Recurrence matrix, redefine/use as necessary. '''
    m = np.zeros((len(x)+1, len(y)+1))
    ''' Recurrence matrix, redefine/use as necessary. '''
    i_x = np.zeros((len(x)+1, len(y)+1))
    ''' Recurrence matrix, redefine/use as necessary. '''
    i_y = np.zeros((len(x)+1, len(y)+1))
    ''' Traceback matrix, redefine/use as necessary. '''
    tb = np.zeros((len(x)+1, len(y)+1, 3))
    tb[:, :, :] = -1

    neg_inf = float('-inf')

    m[0][0] = 0
    i_x[0][0] = 0
    i_y[0][0] = 0
    
    MID = 0
    BOT = 1
    TOP = 2

    for i in range(1, len(x)+1):
        m[i][0] = neg_inf
        i_x[i][0] = i_x[i-1, 0] - e
        i_y[i][0] = neg_inf
        tb[i][0][BOT] = BOT

    for j in range(1, len(y)+1):
        m[0][j] = neg_inf 
        i_x[0][j] = neg_inf
        i_y[0][j] = i_y[0, j-1] - e
        tb[0][j][TOP] = TOP

    # slide 35 in lecture 05
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            from_mid = m[i-1][j-1] + s[x[i-1]][y[j-1]]
            from_bot = i_x[i-1][j-1] +  s[x[i-1]][y[j-1]]
            from_top = i_y[i-1][j-1] + s[x[i-1]][y[j-1]]
            m_score = max(from_mid, from_bot, from_top)
            m_pointer = MID if (m_score == from_mid) else (BOT if (m_score == from_bot) else TOP)
            m[i][j] = m_score

            from_mid = m[i-1][j] - d
            from_bot = i_x[i-1][j] - e
            x_score = max(from_mid, from_bot)
            x_pointer = MID if (x_score == from_mid) else BOT
            i_x[i][j] = x_score

            from_mid = m[i, j-1] - d
            from_top = i_y[i][j-1] - e
            y_score = max(from_mid, from_top)
            y_pointer = MID if (y_score == from_mid) else TOP
            i_y[i][j] = y_score

            tb[i][j][MID] = int(m_pointer)
            tb[i][j][BOT] = int(x_pointer)
            tb[i][j][TOP] = int(y_pointer)

    score = max(m[len(x)][len(y)], i_x[len(x)][len(y)], i_y[len(x)][len(y)])
    tb[len(x)][len(y)][MID] = tb[len(x)][len(y)][TOP] = tb[len(x)][len(y)][BOT] = MID if (score == m[len(x)][len(y)]) else (BOT if (score == i_x[len(x)][len(y)]) else TOP)
    a_x, a_y = traceback(x, y, tb)
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
        description='Calculate sequence alignments for two sequences with an affine gap penalty.')
    parser.add_argument('-f', action="store", dest="f", type=str, required=True)
    parser.add_argument('-s', action="store", dest="s", type=str, required=True)
    parser.add_argument('-d', action="store", dest="d", type=float, required=True)
    parser.add_argument('-e', action="store", dest="e", type=float, required=True)

    args = parser.parse_args()
    fasta_file = args.f
    score_matrix_file = args.s
    d = args.d
    e = args.e

    with open(fasta_file) as f:
        _, x, _, y = [line.strip() for line in f.readlines()]
    with open(score_matrix_file) as f:
        s = json.loads(f.readlines()[0])

    score, (a_x, a_y) = affine_sequence_alignment(x, y, s, d, e)
    print("Alignment:")
    print_alignment(a_x, a_y)
    print("Score: " + str(score))


if __name__ == "__main__":
    main()
