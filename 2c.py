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
    pass


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
    m = None
    ''' Recurrence matrix, redefine/use as necessary. '''
    i_x = None
    ''' Recurrence matrix, redefine/use as necessary. '''
    i_y = None
    ''' Traceback matrix, redefine/use as necessary. '''
    t = None

    a_x, a_y = traceback(x, y, t)
    ''' Complete this function. '''
    pass


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
