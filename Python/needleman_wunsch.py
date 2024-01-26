#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt

def check_match(a, b, match_score, mismatch_score):
    if (a == b):
        return match_score
    else:
        return mismatch_score

def generate_matrix(seq1, seq2, match_score, mismatch_score, gap_penalty):
    sim_matrix = np.zeros((len(seq2) + 1, len(seq1) + 1), dtype=int)
    for i in range(len(seq2) + 1):
        sim_matrix[i, 0] = gap_penalty * i
    for j in range(len(seq1) + 1):
        sim_matrix[0, j] = gap_penalty * j
    for seq1_col in range(len(seq1)):
        for seq2_row in range(len(seq2)):
            match = sim_matrix[seq2_row, seq1_col] + check_match(seq2[seq2_row], seq1[seq1_col], match_score, mismatch_score)
            delete = sim_matrix[seq2_row, seq1_col + 1] + gap_penalty
            insert = sim_matrix[seq2_row + 1, seq1_col] + gap_penalty
            sim_matrix[seq2_row + 1, seq1_col + 1] = max(match, delete, insert)
    return sim_matrix, sim_matrix[-1,-1]

def pairwise_align(seq1, seq2, match_score, mismatch_score, gap_penalty):
    sim_matrix,score = generate_matrix(seq1, seq2, match_score, mismatch_score, gap_penalty)
    AlignmentA = ''
    AlignmentB = ''
    col_seq1 = len(seq1)
    row_seq2 = len(seq2)

    while (col_seq1 > 0 or row_seq2 > 0):
        if (col_seq1 > 0 and row_seq2 > 0):
            diag = sim_matrix[row_seq2, col_seq1] == sim_matrix[row_seq2 - 1, col_seq1 - 1] + check_match(seq2[row_seq2 - 1], seq1[col_seq1 - 1], match_score, mismatch_score)
            delete = sim_matrix[row_seq2, col_seq1] == sim_matrix[row_seq2 - 1, col_seq1] + gap_penalty
            insert = sim_matrix[row_seq2, col_seq1] == sim_matrix[row_seq2, col_seq1 - 1] + gap_penalty
            if (insert):
                AlignmentA += seq1[col_seq1 - 1]
                AlignmentB += '-'
                col_seq1 -= 1
            elif (delete):
                AlignmentA += '-'
                AlignmentB += seq2[row_seq2 - 1]
                row_seq2 -= 1
            elif (diag):
                AlignmentA += seq1[col_seq1 - 1]
                AlignmentB += seq2[row_seq2 - 1]
                col_seq1 -= 1
                row_seq2 -= 1
            else:
                raise Exception('Stuck in an infinite loop')
        elif (col_seq1 > 0):
            AlignmentA += seq1[col_seq1 - 1]
            AlignmentB += '-'
            col_seq1 -= 1
        elif (row_seq2 > 0):
            AlignmentA += '-'
            AlignmentB += seq2[row_seq2 - 1]
            row_seq2 -= 1
        else:
            raise Exception('Stuck in an infinite loop')
    AlignmentA = AlignmentA[::-1]
    AlignmentB = AlignmentB[::-1]

    return AlignmentA, AlignmentB, score


parser = argparse.ArgumentParser()

parser.add_argument(
    '-s',
    '--sequences',
    help = 'Enter the sequences to be compared',
    nargs = 2
)

parser.add_argument(
    '-match',
    '--match_score',
    help = 'Enter the Score given on a match, default is +1',
    default = 1
)

parser.add_argument(
    '-mismatch',
    '--mismatch_score',
    help = 'Enter the Score given on a mismatch, default is -1',
    default = -1
)

parser.add_argument(
    '-g',
    '--gap_penalty',
    help = 'Enter the score given on an insertion or deletion, default is -2',
    default = -2
)

args = parser.parse_args()

if args.sequences:
    seq1 = args.sequences[0]
    seq2 = args.sequences[1]
    match_score = int(args.match_score)
    mismatch_score = int(args.mismatch_score)
    gap_penalty = int(args.gap_penalty)

    align_seq1, align_seq2, score = pairwise_align(seq1, seq2, match_score, mismatch_score, gap_penalty)

    print('----------Sequences---------------')
    print()
    print('Sequence 1: ', seq1)
    print('Sequence 2: ', seq2)
    print()
    print('----------Alignment---------------')
    print()
    print('Alignment 1: ', align_seq1)
    print('Alignment 2: ', align_seq2)
    print()
    print('Score:', score)



